#ifndef TENSOR_FORWARD_HPP
#define TENSOR_FORWARD_HPP

#include "tensor/config.hpp"
#include "tensor/traits.hpp"

#include <cstddef>
#include <boost/multi_array.hpp>

namespace tensor {
namespace detail {

    template<class A>
    struct const_tensor_base;

    template<class A>
    struct tensor_base;

}
}

namespace tensor {

    template<size_t N, typename T = double, class A = boost::multi_array<T,N> >
    struct Tensor;

    template<class E>
    struct expression;

    template<class A>
    struct tensor_ref;

}

#endif // TENSOR_FORWARD_HPP
