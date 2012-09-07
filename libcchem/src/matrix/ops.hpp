#ifndef _MATRIX_OPS_HPP_
#define _MATRIX_OPS_HPP_

namespace matrix {

    template<class M>
    void assign(M &A, typename M::value_type s) {
	for (uint j = 0; j < A.size2(); ++j) {
	    for (uint i = 0; i <  A.size1(); ++i) {
		A(i,j) = s;
	    }
	}
    }

    template<class M, class E>
    typename disable_if< is_arithmetic<E> >::type 
    assign(M &A, const E &e) {
	for (uint j = 0; j < A.size2(); ++j) {
	    for (uint i = 0; i <  A.size1(); ++i) {
		A(i,j) = e(i,j);
	    }
	}
    }

    template<class Matrix>
    void symmeterize(Matrix &A) {
	for (uint j = 0; j < A.size2(); ++j) {
	    for (uint i = j+1; i < A.size1(); ++i) {
		A(i,j) += A(j,i);
		A(j,i) = A(i,j);
	    }
	}
    }

    template<class Matrix>
    void zero(Matrix &A) {
	for (uint j = 0; j < A.size2(); ++j) {
	    for (uint i = 0; i < A.size1(); ++i) {
		A(i,j) = typename Matrix::value_type(0);
	    }
	}
    }

}

#endif /* _MATRIX_OPS_HPP_ */
