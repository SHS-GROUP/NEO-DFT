#ifndef TENSOR_SWAP_HPP
#define TENSOR_SWAP_HPP

#include <algorithm>
#include "tensor/view.hpp"

namespace tensor {

    template<class A, class B, class R>
    void swap(tensor_view<A,R> a, tensor_view<B,R> b) {
	typename tensor_view<A,R>::iterator a_ = a.begin();
	typename tensor_view<B,R>::iterator b_ = b.begin();
	while (a_ != a.end()) {
	    std::swap(*a_++, *b_++);
	}
    }

}

#endif // TENSOR_SWAP_HPP
