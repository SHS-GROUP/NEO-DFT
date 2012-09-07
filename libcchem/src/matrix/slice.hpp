#ifndef _MATRIX_SLICE_HPP_
#define _MATRIX_SLICE_HPP_

#include <boost/array.hpp>

namespace matrix {

    template<size_t N = 2>
    struct slice {
	typedef int index_type;
	typedef size_t size_type;
	typedef boost::array<size_type,2> size_array;
	typedef boost::array<index_type,2> index_array;
	template<class M> class expression;
	slice(const index_array &from, const index_array &to) {
	    for (uint i = 0; i < N; ++i) {
		index_[i] = from[i];
		size_[i] = to[i] - from[i] + 1;
	    }
	}
	template<class M>
	expression<M> operator()(M &A) const {
	    return expression<M>(A, index_, size_);
	}
    private:
	size_array size_;
	index_array index_;
    };

    template<size_t N> template<class M>
    struct slice<N>::expression {
	typedef int index_type;
	typedef typename M::size_type size_type;
	typedef typename M::value_type value_type;
	typedef typename M::reference reference;
	typedef typename M::const_reference const_reference;
	typedef boost::array<size_type,2> size_array;
	typedef boost::array<index_type,2> index_array;
 	expression(M &A) : A_(A) {}
	expression(M &A, const index_array &index, const size_array &size)
	    : A_(A), index_(index), size_(size) {} 
	size_type size1() const { return size_[0]; }
	size_type size2() const { return size_[1]; }
	reference operator()(index_type i, index_type j) {
	    //printf("%i %i\n", i-index_[0], j-index_[1]);
	    return A_(i+index_[0], j+index_[1]);
	}
	const_reference operator()(index_type i, index_type j) const {
	    //printf("%i %i\n", i-index_[0], j-index_[1]);
	    return A_(i+index_[0], j+index_[1]);
	}
    private:
	M &A_;
	index_array index_;
	size_array size_;
    };



}

#endif /* _MATRIX_SLICE_HPP_ */
