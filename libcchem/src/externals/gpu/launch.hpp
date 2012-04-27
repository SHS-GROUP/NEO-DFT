#ifndef GPU_LAUNCH_HPP
#define GPU_LAUNCH_HPP

#include <cstdlib>
#include <cuda.h>

namespace gpu {

    typedef dim3 dim3;

    namespace detail {

	struct launch {
	    explicit launch(dim3 grid = 1, dim3 block = 1, size_t shared = 0)
		: grid(grid), block(block), shared(shared) {}
	    launch operator()(dim3 grid = 1, dim3 block = 1, size_t shared = 0) const {
		return launch(grid, block, shared);
	    }
	    template<class F>
	    void operator[](F f) const;
	private:
	    dim3 grid;
	    dim3 block;
	    size_t shared;
	};

    }
    

    const detail::launch launch = detail::launch();
    
}

#endif //  GPU_LAUNCH_HPP
