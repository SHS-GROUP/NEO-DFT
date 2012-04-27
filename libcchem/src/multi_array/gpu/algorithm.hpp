#ifndef MULTI_ARRAY_GPU_ALGORITHM_HPP
#define MULTI_ARRAY_GPU_ALGORITHM_HPP

#include "foreach.hpp"

#include <cuda.h>
#include <cuda_runtime.h>
#include "multi_array/detail.hpp"
#include "multi_array/algorithm.hpp"
#include <boost/typeof/typeof.hpp>

#include <vector>

#include "boost/cuda/allocator.hpp"
#include "boost/cuda/vector.hpp"

namespace multi_array {
namespace gpu {
namespace algorithm {
	
    struct gpu_tag;

    namespace device {
	    
	template<class F>
	struct kernel {
	    kernel(F f) : stream_(0), f(f) {}
	    kernel(boost::cuda:: stream &stream, F f)
		: stream_(stream), f(f) {}
	    void operator()();
	private:
	    boost::cuda::stream stream_;
	    F f;
	};	

	struct grid;
	struct block;

#ifdef __CUDACC__

	struct grid {
	    __device__
	    grid(int x = blockIdx.x, int y = blockIdx.y,
		 int xsize = gridDim.x, int ysize = gridDim.y)
		: x(x), y(y), xsize(xsize), ysize(ysize) {}
	    const int x, y;
	    const int xsize, ysize;
	};

	struct block 
	{
	    __device__
	    block(int thread = threadIdx.x, int size = blockDim.x)
		: thread(thread), size(size) {}
	    const int thread;
	    const int size;
	};
	
	
	template<class F>
	__global__
	static void global__(F f) {
	    f(grid(), block());
	}

	template<class F>
	void kernel<F>::operator()() {
	    dim3 grid(f.size(0), f.size(1));
	    global__<<<grid, 64, 0, stream_.data<cudaStream_t>()>>>(f);
	    //global__<<<1,1>>>(f);
	}

#endif	    

    }
	
    template<int Order, typename T, typename U, size_t N>
    struct copy_traits {
	typedef multi_array::algorithm::detail::kernel2<
	    Order,
	    detail::multi_array_ref<T,N>,
	    detail::multi_array_ref<U,N>,
	    multi_array::algorithm::detail::copy<gpu_tag> > functor;
	typedef device::kernel<functor> type;
    };
    
    template<int Order, typename T, typename U, size_t N>
    typename copy_traits<Order, T, U, N>::type
    copy(detail::multi_array_ref<T,N> t0,
	 detail::multi_array_ref<U,N> t1) {
	return typename copy_traits<Order, T, U, N>::functor(t0, t1);
    }

    struct pack {
	template<typename T, typename U> struct functor;
	template<typename T, typename U>
	struct kernel {
	    typedef device::kernel<functor<T,U> > type;
	};
	template<class R, typename T, typename U>
	void operator()(const std::vector<R> &range, size_t m,
			 const T *input, size_t ldt,
			 U *output, size_t ldu) const {
	    std::vector<int> index;
	    foreach (const R &r, range) {
		for (BOOST_TYPEOF(r.first) i = 0; i < (r.second - r.first); ++i) {
		    index.push_back(r.first + i);
		}
	    }
	    if (index.empty()) return;
	    index_.assign(index);
	    typedef pack::functor<T,U> F;
	    F f(index_.size(), index_.begin(), m, input, ldt, output, ldu);
	    (device::kernel<F>(stream_, f))();
	}
    private:
	mutable boost::cuda::vector<int> index_;
	mutable boost::cuda::stream stream_;
    public:
	template<typename T, typename U = T>
	struct functor {
	    functor(size_t size, boost::cuda::device_ptr<int> index, size_t m,
		   const T *input, size_t ldt,
		   U *output, size_t ldu)
		: size_(size), index_(index.get()), m_(m),
		  input_(input), ldt_(ldt),
		  output_(output), ldu_(ldu) {}
	    size_t size(int i) const { return ((i == 0) ? size_ : 1); }
	    BOOST_GPU_ENABLED
	    void operator()(const device::grid &grid,
			    const device::block &block)
#ifdef __CUDACC__
	    {
		const T *input = input_ + index_[grid.x]*ldt_;
		U *output = output_ + grid.x*ldu_;
		for (int i = block.thread; i < m_; i += block.size) {
		    output[i] = input[i];
		}
	    }
#else
		;
#endif
	private:
	    size_t size_;
	    int *index_; int m_;
	    const T *input_; int ldt_;
	    U *output_; int ldu_;
	};
    };

} // namespace algorithm
} // namespace gpu

namespace algorithm {
namespace detail {

    template<class T, class U>
    struct functor<copy<gpu::algorithm::gpu_tag>, T, U> {
	template<class I, class O, class B>
	BOOST_GPU_ENABLED
	void operator()(I begin, I end, O output, const B &block) const
#ifdef __CUDACC__
	{
	    for (int i = block.thread; i < (end-begin); i += block.size) {
	    	output[i] = begin[i];
	    }	    
	}
#else
	;
#endif    
    };

}
}
}


#endif // MULTI_ARRAY_GPU_ALGORITHM_HPP
