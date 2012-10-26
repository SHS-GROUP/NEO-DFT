#ifndef ARRAY_ASSIGN_HPP 
#define ARRAY_ASSIGN_HPP

#include <boost/mpl/vector_c.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/utility/swap.hpp>
#include <boost/typeof/typeof.hpp>

#include "array/forward.hpp"
#include "array/detail/tile.hpp"
#include "array/detail/grid.hpp"

#ifdef HAVE_CUDA
#include "cuda.hpp"
#endif // HAVE_CUDA

namespace array {

    template<typename T, typename U, size_t N>
    void assign(adapter<T,N> &a, const adapter<U,N> &b) {
	std::copy(b.data(), b.data() + b.num_elements(), a.data());
    }

#ifdef HAVE_CUDA
    template<typename T, size_t N>
    void assign(adapter<T,N> &a, const adapter<T,N,device_tag> &b) {
	cudaMemcpy(a.data(), b.data(), b.num_elements()*sizeof(T),
		   cudaMemcpyDeviceToHost);
    }
#endif // HAVE_CUDA

}

#endif // ARRAY_ASSIGN_HPP
