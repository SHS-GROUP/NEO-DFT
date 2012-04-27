#ifndef _DISTRIBUTED_OPS_HPP_
#define _DISTRIBUTED_OPS_HPP_

#include "matrix/meta.hpp"
#include "util/unpack.hpp"
#include "distributed.hpp"

namespace distributed {

    using ::matrix::meta_matrix;

    template<class Array, class E>
    void assign(Array &D, const meta_matrix<E> &e) {
	typedef typename Array::index_array index_array;
	typedef typename meta_matrix<E>::index_iterator::indexed indexed;
	foreach (indexed &ab, e.matrix_index()) {
	    typedef typename E::value_type* pointer;
	    UNPACK((int a, b), ab);
	    const E &Mab = e.m(a,b);
	    index_array from = e.matrix_index(a,b);
	    index_array to = from + e.matrix_size(a,b) - 1;
	    index_array ld = e.matrix_size(a,b);
	    pointer src = const_cast<pointer>(Mab.data().begin());

	    //if (from[1] == 16) continue;
	    //std::cout << from << to << ld << "\n";
	    D.put(from, to, src, ld);
	}
    }

    template<class T, size_t N>
    void assign(distributed::array<T,N> &D, T s) {
	D.fill(&s);
    }

}

namespace matrix {

    template<typename T, size_t N, class M>
    void assign(meta_matrix<M> &A, const distributed::array<T,N> &D) {
	typedef typename meta_matrix<M>::index_array index_array;
	typedef typename meta_matrix<M>::index_iterator::indexed indexed;
	foreach (indexed &ab, A.matrix_index()) {
	    typedef typename M::value_type* pointer;
	    UNPACK((int a, b), ab);
	    const M &Mab = A.m(a,b);
	    index_array from = A.matrix_index(a,b);
	    index_array to = from + A.matrix_size(a,b) - 1;
	    index_array ld = A.matrix_size(a,b);
	    pointer src = const_cast< pointer>(Mab.data().begin());

	    if (from[1] == 16) continue;
	    //std::cout << from << to << ld << "\n";
	    //D.get(from, to, src, ld);
	}
    }

}

#endif /* _DISTRIBUTED_OPS_HPP_ */
