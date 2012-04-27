#ifndef GPU_COPY_HPP
#define GPU_COPY_HPP

#include "gpu/forward.hpp"
#include "gpu/runtime.hpp"
#include <iterator>

namespace gpu {

    namespace detail {

	template<class F>
	struct copy;

	template<cudaMemcpyKind K>
	struct copy<copy_tag<K> > {
	    template<class I, class O>
	    copy(I begin, I end, O output) {
		size_t size = sizeof(typename std:: iterator_traits<I>::value_type);
		size *= std::distance(begin, end);
		// std::cout << output << " " << begin << " " << size<<std::endl;
		throw_(cudaMemcpy(output, begin, size, K));
	    }
	};
    

    }


    template<class F, class I, class O>
    void copy(I begin, I end, O output) {
	detail::copy<F>(begin, end, output);
    }

}

#endif // GPU_COPY_HPP
