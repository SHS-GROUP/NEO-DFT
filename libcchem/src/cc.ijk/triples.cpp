#include "config.hpp"
#include "cc/cc.hpp"

#include "cc/triples/triples.hpp"

#include "cc/triples/host.hpp"
#define CC_TRIPLES_THREAD_T3 Host::T3
#include "cc/triples/thread.ijk.hpp"
#undef CC_TRIPLES_THREAD_T3

#ifdef HAVE_CUDA
#include "cc/triples/device.hpp"
#endif

#include <utility>
#include <boost/typeof/typeof.hpp>
#include <boost/numeric/ublas/adaptor.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/multi_array/multi_array_ref.hpp>
#include <boost/array.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include <boost/thread/thread.hpp>
#include <boost/progress.hpp>

#include "foreach.hpp"
#include "array/h5.hpp"

#include "blas.hpp"
#include "omp.hpp"
#include "runtime.hpp"
#include "utility/timer.hpp"

#include "boost/utility/profiler.hpp" 

namespace cc {
namespace triples {
namespace detail {
    
    struct Threads {

	explicit Threads(const runtime &rt)
	    : rt_(rt) {}

	template<class V>
	Result
	operator()(size_t no, size_t nv,
		   const Array<double> &t1,
		   const Array<double> &t2_abij,
		   const Array<double> &t2_aibj,
		   const Array<double> &vvvo,
		   const Array<double> &vvoo,
		   const Array<double> &ovoo,
		   const V &eh, const V &ep,
		   double *work) {

	    Data data(no, nv);

	    data["t1(ai)"] = &t1;
	    data["t2(abij)"] = &t2_abij;
	    data["t2(aibj)"] = &t2_aibj;
	    data["V(abci)"] = &vvvo;
	    data["V(abij)"] = &vvoo;
	    data["V(iajk)"] = &ovoo;

	    data.eh.resize(no);
	    data.ep.resize(nv);
	    std::copy(eh.begin(), eh.end(), data.eh.begin());
	    std::copy(ep.begin(), ep.end(), data.ep.begin());

	    Run run(no, nv);

	    std::vector<int> devices;
#ifdef HAVE_CUDA
	    devices = rt_.get< std::vector<int> >("gpu::devices");
	    boost::ptr_vector<Device::Thread::Team> device;
	    for (size_t i = 0; i < devices.size(); ++i) {
		std::cout << "device: " << i << " " << devices[i] << std::endl;
		device.push_back(new Device::Thread::Team(1, data.no, data.nv));
		run.thread(boost::type<Device>(),
			   Device::Thread(&device[i], i, 0), data);
	    }
#endif // HAVE_CUDA

	    boost::ptr_vector<Host::Thread::Team> teams;
	    int tc = rt_.get<int>("team::count");
	    int ts = rt_.get<int>("team::size");
	    ts -= (devices.size() + tc-1)/tc;
	    ts = std::max(ts,0);
	    for (size_t t = 0; t < tc; ++t) {
		//size_t id = 1;
		size_t dims[] = { nv, nv, nv };
		Host::T3 t3(((t == 0) ? work : 0), dims);
		teams.push_back(new Host::Thread::Team(ts, data.no, data.nv, t3));
		for (size_t i = 0; i < ts; ++i) {
		    std::cout << "team/thread: " 
			      << t << "/" << i << std::endl;
		    run.thread(boost::type<Host>(),
			       Host::Thread(&teams[t], i), data);
		}
	    }
	    //run(boost::type<Host>(), Host::Thread(&host, 0), data);

	    run.join_all();
	    return run.result;
	}

    private:
	runtime rt_;

	struct Run : boost::thread_group, boost::noncopyable {
	    Result result;
	private:
	    boost::mutex mutex_;
	    Task task_;
	    boost::progress_display progress_;

	public:
	    Run(size_t no, size_t nv)
		: task_(no), progress_(task_.size()) {
		result.u1.resize(nv,no);
		std::cout << "tasks: " << nv*task_.size()*4 << std::endl;
	    }

	    template<class T>
	    void thread(boost::type<T>, const typename T::Thread &t,
			const Data &data) {
		this->add_thread(new boost::thread
				 (boost::ref(*this), boost::type<T>(),
				  (t), boost::ref(data)));		
	    }	    

	    template<class T>
	    void operator()(boost::type<T>, const typename T::Thread &t,
			    const Data &data) {
		omp::set_num_threads(1);

		T thread(t);
		BOOST_AUTO(&team, t.team());
		bool must_free = false;
		if (thread.master && !team.t3.data()) {
		    std::cout << "allocate t3" << std::endl;
		    size_t nv = data.nv;
		    size_t shape[] = { nv, nv, nv };
		    BOOST_AUTO(t3, T::template make_array<3>(shape));
		    team.t3.swap(t3);
		    must_free = true;
		}
		team.wait();
		
		while (true) {
		    team.wait();
		    if (thread.master) {
			team.task = task_.next();
		    }
		    team.wait();
		    if (!team.task) break;

		    thread(data, team.t3, team.task);
		    if (thread.master) {
			boost::lock_guard<boost::mutex> lock(mutex_);
			++progress_;
		    }
		}
		team.wait();
		if (must_free) T::free(team.t3.data());
		{
		    boost::lock_guard<boost::mutex> lock(mutex_);
		    result.C += thread.C;
		    result.u1 += thread.u1;
		    //std::cout << thread.id << ":  " << thread.C.etd << std::endl;
		}
	    }
	};

    };

} // namespace triples
} // namespace detail
} // namespace cc


cc::Triples::Correction
cc::Triples::operator()(size_t no, size_t nv,
			const Data &data,
			const vector &eh, const vector &ep,
			double* t3) {

    struct {
	const Array<double> & operator()(const std::string &k) const {
	    return data.find(k)->second;
	}
	const Data &data;
    } find = { data };
    

    utility::timer t;

    size_t dims[] = { nv, no, nv, no };
    Array<double, 0, H5::H5File> t2_aibj(4, dims, "t(a,i,b,j)");
    {
	assert(no <= nv);
	triples::detail::Matrix vo(nv, no);
	boost::multi_array_ref<double,3> t2_(t3, boost::extents[no][nv][nv]);
	for (int j = 0; j < int(no); ++j) {
	    size_t r1[] = { 0, 0, 0, j };
	    size_t r2[] = { nv, nv, no, j+1 };
	    find("t2").get(t2_.origin(), r1, r2);
	    for (int b = 0; b < int(nv); ++b) {
		for (int i = 0; i < int(no); ++i) {
		    namespace ublas = boost::numeric::ublas;
		    ublas::column(vo,i) = ublas::make_vector(nv, &t2_[i][b][0]);
		}
		size_t s1[] = { 0, 0, b, j };
		size_t s2[] = { nv, no, b+1, j+1 };
		t2_aibj.put(vo.data().begin(), s1, s2);
	    }
	}
    }

    boost::utility::global_profiler().clear();
    //boost::multi_array_ref<double,3> t3_(t3, boost::extents[nv][nv][nv]);

    runtime rt("CC");
    triples::detail::Threads threads(rt);
    //(runtime::num_threads("CC"), runtime::gpu_devices());
    const triples::detail::Result &result =
	threads(no, nv,
		find("t1"), find("t2"), t2_aibj,
		find("vvvo"), find("vvoo"), find("ovoo"),
		eh, ep, t3);

    Correction C = result.C;
    triples::detail::Matrix t1(nv, no, 0);
    const size_t r1[] = { 0, 0 };
    const size_t r2[] = { nv, no };
    find("t1").get(t1.data().begin(), r1, r2);
    double ets = 0;
    for (int i = 0; i < int(nv*no); ++i) {
	ets += t1.data()[i]*result.u1.data()[i];
    }
    C.ets += 2*ets;

    std::cout <<  "CC triples: " << t << std::endl;

    std::cout << boost::utility::global_profiler() << std::endl;

    return C;

}
