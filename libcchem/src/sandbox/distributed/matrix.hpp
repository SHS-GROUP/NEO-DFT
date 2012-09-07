#ifndef _DISTRIBUTED_MATRIX_HPP_
#define _DISTRIBUTED_MATRIX_HPP_

#include "distributed/array.hpp"

namespace distributed {

    template< typename T>
    struct matrix : public  distributed::array<T,2> {
	typedef distributed::array<T,2> base;
	typedef typename base::value_type value_type;
	typedef typename base::index_type index_type;
	typedef typename base::size_type size_type;
	typedef typename base::index_array index_array;

	typedef typename base::const_pointer const_pointer;

	typedef typename base::element element;
	typedef typename base::const_element const_element;

	matrix(size_t m, size_t n) : base(index_array(m,n)) {}

	template<class E>
	void operator=(const E &e) { base::operator=(e); }

	size_type size1() const  { return this->dims()[0]; }

	size_type size2() const { return this->dims()[1]; }

	element operator()(index_type i, index_type j) {
	    return element(*this, index_array(i,j));
	}

	const_element operator()(index_type i, index_type j) const {
	    return const_element(*this, index_array(i,j));
	}

    };

}

#endif /* _DISTRIBUTED_MATRIX_HPP_ */
