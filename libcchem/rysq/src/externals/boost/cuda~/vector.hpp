#ifndef BOOST_CUDA_VECTOR2_HPP
#define BOOST_CUDA_VECTOR2_HPP

#include <vector>
#include "boost/cuda/allocator.hpp"
#include "boost/cuda/array.hpp"

namespace boost {
namespace cuda {

    template<typename T>
    struct mapped_vector {
	typedef mapped_allocator<T> allocator;
	typedef std::vector<T, allocator> type;
    };

}
}

#endif // BOOST_CUDA_VECTOR2_HPP
