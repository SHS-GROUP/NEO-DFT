#ifndef _MATRIX_PERMUTE_HPP_
#define _MATRIX_PERMUTE_HPP_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>

namespace matrix {

    template<typename T>
    struct permutation {
	typedef permutation self;

	template<class M, class Permutation = self>
	struct expression { // : ublas::matrix_expression< matrix<Matrix, Permutation> > {
	    typedef expression self;
	    typedef typename M::size_type size_type;
	    typedef typename M::reference reference;
	    typedef typename M::const_reference const_reference;
	    expression(M &A, const Permutation &P) : A(A), P(P) {}
	    size_type size1() const { return A.size1(); }
	    size_type size2() const { return A.size2(); }
	    reference operator()(size_type i, size_type j) {
		return A(P[i], P[j]);
	    }
	    const_reference operator()(size_type i, size_type j) const {
		return A(P[i], P[j]);
	    }
	    template<class E>
	    self& operator=(const E &e) { assign(*this, e); return *this; }
	private:
	    M &A;
	    const Permutation &P;
	};

	const std::vector<T> &index;
	permutation(const std::vector<T> &index) : index(index) {}
	index_type operator[](index_type i) const { return index[i]; }

	template<class M>
	expression<M,self> operator()(M &A) const {
	    return expression<M,self>(A, *this);
	}
    };



}

#endif /* _MATRIX_PERMUTE_HPP_ */
