#ifndef HF_DEVICE_HPP
#define HF_DEVICE_HPP

#include <rysq.hpp>

#if RYSQ_CUDA
#define HF_CUDA 1

#include "hf/hf.hpp"
#include "hf/thread.hpp"
#include "hf/histogram.hpp"
#include "adapter/rysq.hpp"
#include "omp.hpp"
#include "cuda.hpp"

#include <memory>
#include <boost/ptr_container/ptr_map.hpp>


#include "boost/utility/profiler.hpp"

namespace hf {

    template<class Matrix, class Quartets>
    struct Thread<Matrix, Quartets>::Device {
    private:

	Thread<Matrix, Quartets> host_;

	struct Stream : cuda::stream {
	    rysq::cuda::Quartets quartets;
	    ::cuda::event event;
	    explicit Stream(cudaStream_t data = NULL)
		: cuda::stream(data) {}
	    operator rysq::cuda::Stream() {
		return rysq::cuda::Stream((cudaStream_t)*this);
	    }
	};

	typedef boost::multi_array<rysq::cuda::density_matrix,2> density_array_type;
	typedef boost::multi_array<rysq::cuda::fock_matrix,2> fock_array_type;

    public:

	template<class Blocks, class Centers, class Density, class Fock,
		 class S>
	void operator()(const Blocks &blocks, const Centers &centers,
			const Density &D, lockable<Fock> &F,
			const S &Screen, work_queue &queue,
			const Parameters &parameters) {

	    BOOST_PROFILE_LINE;
	    timer t;
	    Histogram histogram;

	    boost::array<size_t,2> shape = {{ D.array()[0].size(),
					      D.array()[1].size() }};
	    density_array_type density_array(shape);
	    fock_array_type fock_array(shape);

	    for (size_t j = 0; j < shape[1]; ++j) {
		for (size_t i = 0; i < shape[0]; ++i) {
		    size_t size1 = D.m(i,j).size1();
		    size_t size2 = D.m(i,j).size2();
		    size_t block1 = D.m(i,j).block1();
		    size_t block2 = D.m(i,j).block2();
		    boost::array<int,2> ij = {{ i, j}};

		    fock_array(ij).resize(size1, size2, block1, block2);
		    density_array(ij).resize(size1, size2, block1, block2);

		    fock_array(ij).clear();
		    density_array(ij).assign(D.m(i,j));
		}
	    }

	    std::vector<rysq::Int4> quartets;

	    using  rysq::hf::matrix_ptr_set;

	    struct {
		rysq::cuda::Centers centers;
		
		std::auto_ptr<matrix_ptr_set<rysq::cuda::fock_matrix> > F;
		std::auto_ptr<matrix_ptr_set<rysq::cuda::density_matrix> > D;

		rysq::cuda::Mutex mutex;
		std::vector<Stream> streams;
		Stream& stream() {
		    foreach (Stream& s, streams) {
			if (s.query()) return s;
		    }
		    return streams.at(next_++%streams.size());
		}
	    private:
		size_t next_;
	    } device;

	    {
		size_t num_shells = 0;
		typedef BOOST_TYPEOF(blocks[0]) B;
		foreach (const B &b, blocks) {
		    num_shells += b.size();
		}
		//std::cout << num_shells << std::endl;
		device.mutex.resize(num_shells);
		device.centers.assign(centers);
		device.streams.resize(8);
		foreach (Stream &s, device.streams) {
		    s.create();
		    s.event.create();
		}
	    }

	    boost::ptr_map<
		Basis::Shell::Data::key_type, rysq::cuda::Shell> cache;

	    rysq::cuda::Context context;

	    BOOST_PROFILE_LINE;
	    // Stream stream;

	    while (true) {

		BOOST_PROFILE_LINE;

		boost::array<int,4> task;
		{
		    // BOOST_PROFILE_LINE;
		    try {
			work_queue::value_type t = queue.pop();
			for (int i = 0; i < 4; ++i) {
			    task[i] = t[i];
			}
		    }
		    catch(std::exception&) { break; }
		}
		std::swap(task[1], task[2]); // ikjl -> ijkl
		BOOST_AUTO(screen, Screen(blocks, task, (1<<15) ));

		struct {
		    int base[4], block[4];
		} data;

		boost::array<rysq::cuda::Shell*,4> Q = {{}};

		for (int i = 0; i < 4; ++i) {
		    const Basis::Block &b = blocks[task[i]];
		    BOOST_AUTO(k, b.shell().key());
		    if (!cache.count(k)) {
			typedef rysq::cuda::Shell Shell;
			Shell *s = new Shell(adapter::rysq::Shell(b.shell()));
			cache.insert(k, s);
		    }
		    Q[i] = &cache.at(k);
		}

		Transpose transpose;
		{
		    int value = 0;
		    if (Q[0]->size() < Q[1]->size())
			value = value | Transpose::BRA;
		    if (Q[2]->size() < Q[3]->size())
			value = value | Transpose::KET;
		    if (Q[0]->L() + Q[1]->L() < Q[2]->L() + Q[3]->L())
			value = value | Transpose::BRAKET;
		    transpose = Transpose(value);
		    //std::cout << value << std::endl;
		}

		task = transpose(task);
		Q = transpose(Q);

		int size = 1;
		for (int i = 0; i < 4; ++i) {
		    const Basis::Block &b = blocks[task[i]];
		    size *= Q[i]->size();
		    BOOST_AUTO(it, blocks.begin());
		    data.base[i] = b.start();
		    data.block[i] = std::distance(it, it + task[i]);
		}

		{
		    BOOST_PROFILE_LINE;
		    device.D.reset
			(new matrix_ptr_set<rysq::cuda::density_matrix>(data.base));
		    device.F.reset
			(new matrix_ptr_set<rysq::cuda::fock_matrix>(data.base));
		}

		rysq::cuda::Fock f(*Q[0], *Q[1], *Q[2], *Q[3], context);

		foreach (const size_t (&index)[2], rysq::index_list()) {
		    size_t i = index[0], j = index[1];
		    size_t im = data.block[i];
		    size_t jm = data.block[j];
		    // std::cout << im << " " << jm << std::endl;
		    device.F->get(i,j) = &fock_array[im][jm];
		    device.D->get(i,j) = &density_array[im][jm];
		}

		BOOST_PROFILE_LINE;

		while (true) {
		    Stream &stream = device.stream();
		    BOOST_PROFILE_LINE;
		    if (f) {// && size == 108) {
			BOOST_PROFILE_LINE;
			{
			    BOOST_PROFILE_LINE;
			    stream.event.synchronize();
			    screen(quartets, transpose);
			    if (quartets.empty()) break;
			}
			{
			    BOOST_PROFILE_LINE;
			    stream.synchronize();
			    stream.quartets.assign(quartets, stream);
			    stream.event.record(stream);
			}
			{
			    //stream.synchronize();
			    BOOST_PROFILE_LINE;
			    timer t;
			    f(device.centers, stream.quartets,
			      *device.D, *device.F, device.mutex,
			      rysq::cuda::Fock::Parameters
			      (parameters.cscale, parameters.xscale,
			       parameters.cutoff),
			      stream);
			    //stream.synchronize();
			    histogram.gpu[size] += t;
			}
		    }
		    else {
		    	// BOOST_PROFILE_LINE;
			timer t;
			screen(quartets, transpose);
 		    	host_(blocks, centers, D, F,
			      task, quartets, parameters);
			histogram.cpu[size] += t;
		    }
		    //BOOST_PROFILE_LINE;
		    if (quartets.empty()) break;
		    quartets.clear();
		}


	    }	    

	    BOOST_PROFILE_LINE;

	    ::cuda::synchronize();

	    for (size_t j = 0; j < fock_array[1].size(); ++j) {
	    	for (size_t i = 0; i < fock_array[0].size(); ++i) {
		    BOOST_AUTO(const &f, rysq::cuda::host(fock_array[i][j]));
		    omp::scoped_lock lock(F.lock(i,j));
	    	    F.data().m(i,j) += f;
	    	}
	    }

	    //rysq::cuda::finish();

	    // std::cout << "Thread histogram:" << std::endl;
	    // std::cout << histogram << std::endl;
	    // std::cout << "total: " << histogram.total() << std::endl;

	}


    };

} // namespace hf

#endif /* RYSQ_CUDA */

#endif /* HF_DEVICE_HPP */
