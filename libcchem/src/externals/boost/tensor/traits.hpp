#ifndef TENSOR_TRAITS_HPP
#define TENSOR_TRAITS_HPP

#include "tensor/storage/storage.hpp"
#include "tensor/storage/multi_array.hpp"

namespace tensor {
namespace detail {

    template<class AE>
    struct traits : storage::array<AE> {
	typedef typename storage::array<AE>::reference result_type;
	typedef typename storage::array<AE>::const_reference const_result_type;
    };

    

}
}

#endif // TENSOR_TRAITS_HPP
