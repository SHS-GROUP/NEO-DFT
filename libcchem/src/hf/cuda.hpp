#ifndef _HF_CUDA_HPP_
#define _HF_CUDA_HPP_

#include "hf/hf.hpp"
#include "hf/thread.hpp"

#include <memory>
#include <boost/thread.hpp>
#include <boost/thread/exceptions.hpp>
#include <boost/multi_array.hpp>
#include <boost/logic/tribool.hpp>

#include "adapter/rysq.hpp"
#include <rysq.hpp>

#if RYSQ_CUDA
#define HF_CUDA 1

namespace hf {
namespace detail {
	
    struct gpu_fock {
    private:
	struct thread_data {
	    thread_data()
		: state(boost::logic::indeterminate) {}
	    void wait(boost::unique_lock<boost::mutex> &lock) {
		condition_.wait(lock);
	    }
	    void notify() { condition_.notify_one(); }
	    //boost::barrier barrier;
	    boost::logic::tribool state; 
	    bool restart;
	    boost::mutex mutex; 
	    rysq::cuda::Fock *fock;
	    boost::array<int,4> base, block; 
	    rysq::cuda::Fock::Quartets *quartets;
	private:
	    boost::condition_variable condition_;
	};


    public:

	typedef boost::multi_array<rysq::cuda::density_matrix,2> density_array_type;
	typedef boost::multi_array<rysq::cuda::fock_matrix,2> fock_array_type;

	template<class B, class C, class S, class D, class F, class Q>
	void operator()(const B &blocks, const C &centers, const S &Screen,
			const D &density, lockable<F> &fock, const Parameters &p,
			Q &queue, int device_id) {

	    //timer t; 

	    rysq::cuda::Fock::Quartets screened;
	    rysq::cuda::Fock::Parameters
		parameters(p.cscale, p.xscale, Screen.cutoff());

	    rysq::cuda::initialize(device_id);

	    thread_data data;
	    using boost::ref;
	    boost::thread device(*this, ref(centers), ref (density), ref(fock),
				 parameters, device_id, ref(data));

	    while (true) {
		typename Q::value_type task;
		try { task = queue.pop(); }
		catch(std::exception&) { break; }

		BOOST_AUTO(ib, &blocks[task[0]]);
		BOOST_AUTO(kb, &blocks[task[1]]);
		BOOST_AUTO(jb, &blocks[task[2]]);
		BOOST_AUTO(lb, &blocks[task[3]]);
		const BOOST_TYPEOF(ib) block[] = { ib, jb, kb, lb };

		//if (!screen.test(block)) return true;

		typedef adapter::rysq::Shell Shell;
		Shell a(block[0]->shell());
		Shell b(block[1]->shell());
		Shell c(block[2]->shell());
		Shell d(block[3]->shell());

		rysq::cuda::Fock f(rysq::Quartet<rysq::Shell>(a,b,c,d));
		if (!f) {
		    queue.backlog().push(task);
		    continue;
		}
		BOOST_AUTO(screen, Screen(block, (1<<21)));

		screened.clear();
		screen(screened);

		data.fock = &f;
		for (int i = 0; i < 4; ++i) {
		    data.base[i] = block[i]->start();
		    data.block[i] = std::distance(&blocks[0], block[i]);
		}
		data.restart = false;
		while (!screened.empty()) {
		    {
			boost::lock_guard<boost::mutex> lock(data.mutex);
			data.quartets = &screened;
			data.state = true;
			data.notify();
		    }
		    {
			boost::unique_lock<boost::mutex> lock(data.mutex);
			while (data.state) { data.wait(lock); } 
		    }
 		    screened.clear();
		    screen(screened);
		    data.restart = true;
		}
	    }

	    //std::cout << "GPU thread time: " << t << std::endl;

	    {
		boost::lock_guard<boost::mutex> lock(data.mutex);
		data.state = false;
	    }
	    data.notify();
	    device.join();

	    rysq::cuda::finish();
	}

	template<class C, class D, class F, class P>
	void operator()(const C &centers, const D &density, lockable<F> &fock,
			const P &parameters, int device_id, thread_data &data) {
	    
	    rysq::cuda::thread thread(device_id);

	    boost::array<size_t,2> shape = {{ density.array()[0].size(),
					      density.array()[1].size() }};

	    density_array_type density_array(shape);
	    fock_array_type fock_array(shape);

	    for (size_t j = 0; j < shape[1]; ++j) {
		for (size_t i = 0; i < shape[0]; ++i) {
		    size_t size1 = density.m(i,j).size1();
		    size_t size2 = density.m(i,j).size2();
		    size_t block1 = density.m(i,j).block1();
		    size_t block2 = density.m(i,j).block2();
		    boost::array<int,2> ij = {{ i, j}};

		    fock_array(ij).resize(size1, size2, block1, block2);
		    density_array(ij).resize(size1, size2, block1, block2);

		    fock_array(ij).clear();
		    density_array(ij).assign(density.m(i,j));
		}
	    }

	    struct {
		rysq::cuda::Centers centers;
		rysq::cuda::Fock::Quartets quartets;
	    } device;
	    device.centers.assign(centers);

	    rysq::cuda::Fock f;
	    while (true) {
		BOOST_AUTO(base, data.base);
		BOOST_AUTO(block, data.block);
		
		{
		    using boost::logic::indeterminate;
		    
		    boost::unique_lock<boost::mutex> lock(data.mutex);
		    while (indeterminate(data.state)) { data.wait(lock); } 
		    if (!data.state) break;
		    //std::cout << "GPU thread: quartet" << std::endl;

		    if (!data.restart) f.swap(*data.fock);
		    device.quartets.swap(*data.quartets);
		    base = data.base;
		    block = data.block;
		    data.state = indeterminate;
		}
		data.notify();
		
		rysq::hf::matrix_ptr_set<rysq::cuda::fock_matrix> fock(base.elems);
		rysq::hf::matrix_ptr_set<rysq::cuda::density_matrix> density(base.elems);

		foreach (const size_t (&index)[2], rysq::index_list()) {
		    size_t i = index[0], j = index[1];
		    size_t im = block[i];
		    size_t jm = block[j];
		    // std::cout << im << " " << jm << std::endl;
		    fock.get(i,j) = &fock_array[im][jm];
		    density.get(i,j) = &density_array[im][jm];
		}
		f(device.centers, device.quartets, density, fock, parameters);
	    }	    

	    for (size_t j = 0; j < fock_array[1].size(); ++j) {
	    	for (size_t i = 0; i < fock_array[0].size(); ++i) {
	    	    BOOST_AUTO(const &f, rysq::cuda::host(fock_array[i][j]));
	    	    boost::lock_guard<boost::mutex> lock(fock.lock(i,j));
	    	    fock.data().m(i,j) += f;
	    	}
	    }
		
	}
    }; 

} // namespace detail
}

namespace hf {

    struct cuda {

	static std::vector<int> all_devices() {
	    return rysq::cuda::get_devices();
	}

	struct fock_thread_group : boost::noncopyable {
	    boost::thread_group group_;
	    template<class B, class C, class S, class D, class F, class Q>
	    fock_thread_group(const B &blocks, const C &centers, const S &screen,
			      const D &density, lockable<F> &fock,
			      const Parameters &parameters, Q &queue,
			      const std::vector<int> &devices = cuda::all_devices()) {
		for (size_t i = 0; i < devices.size(); ++i) {
		    group_.add_thread
			(fock_thread::create
			 (blocks, centers, screen, density, fock,
			  parameters, queue, devices[i]));
		}
	    }
	    void join_all() { group_.join_all(); }
	};

	struct fock_thread : boost::noncopyable {

	    template<class B, class C, class S, class D, class F, class Q>
	    static
	    boost::thread* create(const B &blocks, const C &centers,
				  const S &screen,
				  const D &density, lockable<F> &fock,
				  const Parameters &parameters,
				  Q &queue, int device = 0) {
		using boost::ref;
		return new 
		    boost::thread(detail::gpu_fock(),
				  blocks, ref(centers), ref(screen), 
				  ref(density), ref(fock), 
				  parameters, ref(queue), device);
	    }

	    template<class B, class C, class S, class D, class F, class Q>
	    fock_thread(const B &blocks, const C &centers, const S &screen,
			const D &density, lockable<F> &fock,
			const Parameters &parameters, Q &queue, int device = 0)
		: thread_(create(blocks, centers, screen, density, fock,
				 parameters, queue, device))
	    {
	    } 
	    void join() { thread_.get()->join(); }

	private:
	    std::auto_ptr<boost::thread> thread_; 
	};
    };

}


#endif /* RYSQ_CUDA */

#endif /* _HF_CUDA_HPP_ */
