#include "cc/triples/device.hpp"

#define CC_TRIPLES_THREAD_T3 Device::T3
#include "cc/triples/thread.ijk.hpp"
#undef CC_TRIPLES_THREAD_T3

// #include <utility>
// #include <boost/typeof/typeof.hpp>
// #include <boost/numeric/ublas/adaptor.hpp>
// #include <boost/numeric/ublas/io.hpp>
// #include <boost/multi_array/multi_array_ref.hpp>
// #include <boost/array.hpp>

// #include <boost/fusion/include/pair.hpp>
// #include <boost/fusion/include/map.hpp>
// #include <boost/fusion/include/intrinsic.hpp>
// #include <boost/fusion/include/for_each.hpp>

// #include <boost/thread/thread.hpp>

// #include "foreach.hpp"
// #include "array/h5.hpp"
// #include "array/permute.hpp"

// #include "blas.hpp"
// #include <cublas.h>

#define BOOST_PROFILE_ENABLE
#include "boost/utility/profiler.hpp"

namespace cc {
namespace triples {
namespace detail {


    void* Device::allocate_impl(size_t size) {
	void *ptr = 0;
	cudaMalloc(&ptr, size);
	return ptr;
    }

    void Device::free_impl(void* ptr) {
	cudaFree(ptr);
    }

    Device::Device(const Thread &thread, size_t ns)
    	: Thread(thread)
    {
    	done = 0;
    	ns = (ns ? ns : 4);//num_streams(Thread::device));
    	cudaSetDevice(Thread::device);
    	cublas::init();
    	cublasInit();

    	std::cout << "thread: " << Thread::id << ", "
    		  << "device: " << Thread::device << " "
    		  << "streams: " << ns << std::endl;
	 
    	int nv = Thread::nv_;
    	int no = Thread::no_;
    	tmp_.resize(nv, nv*ns, false);
    	t2_.resize(nv, nv, false);
    	Vjk_.resize(no, nv, false);
    	Vkj_.resize(no, nv, false);
	vov_.resize(nv, no*nv, false);
	hvvv_.resize(nv, nv*nv, false);
	vvv_.resize(nv, nv*nv, false);
    	foreach(cublas::matrix<double> &m, v_) {
    	    m.resize(nv, std::max(nv,no)*ns, false);
    	}
    	streams_.resize(ns, 0);
    	foreach(cudaStream_t &s, streams_) {
    	    cudaStreamCreate(&s);
    	}
    	//cudaHostAlloc(&H_, nv*nv*ns*sizeof(double), cudaHostAllocDefault);
    }

    Device::~Device() {
    	std::cout << "done: "<< done << std::endl;
    	foreach(cudaStream_t &s, streams_) {
    	    cudaStreamDestroy(s);
    	}
    }


    void Device::ijk1_impl(const Data &data,
    			   const Matrix &tab, const Matrix &taib,
    			   const Matrix &Vjk, const Matrix &Vkj,
    			   Device::T3 &t3,
    			   int i, int j, int k,
    			   std::pair<bool,int> terms,
    			   Queue &queue) {


    	namespace ublas = boost::numeric::ublas;
    	//Layout layout;
    	int nv = Thread::nv_;
    	int no = Thread::no_;

    	if (terms.first) {
    	    t2_.assign(tab);

		{
		    BOOST_PROFILE("get");
    	    	data["V(abci)"]
    	    	    .get(hvvv_.data().begin(),
    	    		 index<4>(0,0,0,k), index<4>(nv,nv,nv,k+1));
		}
		{
		    BOOST_PROFILE("copy::0");
		    vvv_.assign(hvvv_);
		}

    	}
    	if (terms.second) {
    	    Vjk_.assign(Vjk);
    	    Vkj_.assign(Vkj);
	    BOOST_PROFILE("copy::1");
	    vov_.assign(taib);
    	}

	cudaThreadSynchronize();
	//cudaHostRegister(tmp_.data().begin(), nv*nv*streams_.size(), 0);

    	while (true) {
    	    Queue::Tuple w = queue.advance(streams_.size());
    	    int c = w.begin;
    	    int n = w.size;
    	    if (n < 1) break;
    	    done += n;

    	    double b = 0;
    	    ublas::range rv(0, nv), rvn(0, nv*n);
	    ublas::range ro(0, no), ron(0, no*n);

	    cublas::matrix_adaptor<double> t3_(nv, nv*n, &t3[c][0][0]);

    	    BOOST_AUTO(H, ublas::project(tmp_, rv, rvn));

    	    if (terms.first) {
		// {
		//     BOOST_PROFILE("get");
    	    	// data["V(abci)"]
    	    	//     .get(tmp_.data().begin(),
    	    	// 	 index<4>(0,0,c,k), index<4>(nv,nv,c+n,k+1));
		// }
		// {
		//     BOOST_PROFILE("copy::0");
		//     BOOST_AUTO(vv, cublas::project(v_[0], rv, rvn));
		//     vv = ublas::project(tmp_, rv, rvn);
		// }

    	    	for (int s = 0; s < n; ++s) {
    	    	    cublas::set_stream(streams_.at(s));
    	    	    ublas::range rs(s*nv, (s+1)*nv);
    	    	    BOOST_AUTO(vv, cublas::project(v_[0], rv, rs));
    	    	    BOOST_AUTO(t3_s, cublas::project(t3_, rv, rs));
		    BOOST_PROFILE("gemm::0");
    	    	    cublas::gemm(1, vv, t2_, 0, t3_s);
    	    	}
    	    	//cudaThreadSynchronize();
		BOOST_PROFILE("gemm::1");
		cublas::gemm(1, t2_, cublas::project(v_[0], rv, rvn), 1, t3_);
    	    	b = 1;
    	    }
	    
    	    for (int s = 0; s < n; ++s) {
    		cublas::set_stream(streams_.at(s));

    		ublas::range rvs(s*nv, (s+1)*nv);
    		ublas::range ros(no*(c+s), no*(c+s+1));

    		BOOST_AUTO(vo, cublas::project(vov_, rv, ros));

    		if (terms.second == 2) {
		    BOOST_PROFILE("gemm::2");
    		    BOOST_AUTO(t3_s, cublas::project(t3_, rv, rvs));
    		    cublas::gemm(-1, vo, Vjk_, b, t3_);
    		}

    		if (terms.second == 1) {
		    BOOST_PROFILE("gemm::3");
    		    BOOST_AUTO(t3_b, cublas::project(v_[2], rv, rvs));
    		    cublas::gemm(-1, vo, Vkj_, 0, t3_b);
    		}

    	    }		


    	    // if (terms.second != 1) {
    	    //     cublas::host(H) = cublas::project(v_[1], rv, rvn); 
    	    //     BOOST_AUTO(t3_, ublas::make_matrix(nv, nv*n, t3[c].origin(), layout));
    	    //     ublas::noalias(t3_) += H;
    	    // }
    	    // else {
    	    //     cublas::host(H) =
    	    //      	cublas::project(v_[2], rv, ublas::range(0, nv*n));
    	    //     for (int s = 0; s < n; ++s) {
    	    //     	for (int b = 0; b < nv; ++b) {
    	    //     	    BOOST_AUTO(t3_bc,
    	    //     		       ublas::make_vector(nv, t3[b][c+s].origin()));
    	    //     	    ublas::noalias(t3_bc) += ublas::column(tmp_,b+s*nv);
    	    //     	}
    	    //     }
    	    // }
		

    	} // while

	cudaThreadSynchronize();

    }

}
}
}

