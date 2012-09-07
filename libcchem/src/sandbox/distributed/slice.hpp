#ifndef _DISTRIBUTED_SLICE_HPP_
#define _DISTRIBUTED_SLICE_HPP_

#include <boost/utility.hpp>
#include "util/slice.hpp"

namespace distributed {

    namespace slice {

	template<typename T>
	struct range {
	    typedef T index_type;
	    typedef std::vector<index_type> index_vector;
	    typedef size_t size_type;

	    template<typename U, size_t N>
	    range(const U (&from)[N], const U (&to)[N])
		: from_(to_vector(from)), to_(to_vector(to)) {}
	    
	    index_type* c_from() const { return c_array(from_); }
	    index_type* c_to() const { return c_array(to_); }

	    template<typename U, size_t N>
	    static index_vector to_vector(const U (&index)[N]) {
		index_vector v(N);
		std::copy(index, index + N, v.begin());
		return v;
	    }
	    const index_vector& from() const { return from_; }
	    const index_vector& to() const { return to_; }
	    size_type size(index_type i) const { return to_[i] - from_[i] + 1; }
	private:
	    const index_vector from_, to_;
	    static index_type* c_array(const index_vector &v) {
		return const_cast< index_type*>(&v[0]);
	    }
	};

	template<class Array>
	struct operators {
	    typedef typename Array::value_type value_type;
	    typedef typename Array::pointer pointer;
	    typedef typename Array::const_pointer const_pointer;
	    typedef typename Array::index_type index_type;
	    typedef range<index_type> range_type;

	    operators(Array &A, const range_type &range) : A_(A), range_(range) {}
	    void operator=(value_type t) {
		if (t == value_type(0))  apply(&Array::zeroPatch);
		else apply(&Array::fillPatch, t);
	    }
	    void operator+=(value_type t) { apply(&Array::addConstantPatch, t); }
	    void operator-=(value_type t) { operator+=(-t); }
	    void operator*=(value_type t) { apply(&Array::scalePatch, t); }
	    void operator/=(value_type t) { operator*=(value_type(1)/t); }

	protected:
	    void assign(const_pointer p, const index_type *ld) {
		A_.put(from(), to(), const_cast<pointer>(p),
		       const_cast<index_type*>(ld));
	    }
	    template<class F>
	    void apply(F f) { (A_.*f)(from(), to()); }
	    template<class F>
	    void apply(F f, value_type v) { (A_.*f)(from(), to(), &v); }
	private:
	    Array &A_;
	    const range_type &range_;
	    index_type* from() const { return  range_.c_from(); }
	    index_type* to() const { return  range_.c_to(); }
	};


    }

    template<typename T, size_t M, class Array = distributed::array_base<T> >
    struct array_slice :
	slice::range<typename Array::index_type>,
	slice::operators<Array>
    {

	typedef slice::range<typename Array::index_type> range_base;
	typedef slice::operators<Array> operators_base;

	template<typename U, size_t N>
	array_slice(Array &A, const U (&from)[N], const U (&to)[N])
	    : range_base(from, to), operators_base(A, *this) {}

	template<typename U, size_t N>
	array_slice(Array &A, const boost::array<U,N> &from, const boost::array<U,N> &to)
	    : range_base(from.elems, to.elems), operators_base(A, *this) {}

	using operators_base::operator=; 

    };

    template<typename T, class Array>
    struct array_slice<T, 2, Array> :
	slice::range<typename Array::index_type>,
	slice::operators<Array>
    {
	typedef size_t size_type;
	typedef typename Array::index_type index_type;
	typedef typename Array::const_pointer const_pointer;

	typedef slice::range<index_type> range;
	typedef slice::operators<Array> operators_base;

	template<typename U, size_t N>
	array_slice(Array &A, const U (&from)[N], const U (&to)[N])
	    : range(from, to), operators_base(A, *this), N_(N) {}

	template<typename U, size_t N>
	array_slice(Array &A, const boost::array<U,N> &from, const boost::array<U,N> &to)
	    : range(from.elems, to.elems), operators_base(A, *this), N_(N) {}

	using operators_base::operator=; 

	template<class E>
	void operator=(const ublas::matrix_expression<E> &AE) {
	    assert(size1() == AE().size1());
	    assert(size2() == AE().size2());
	    const_pointer p = ublas::matrix<T>(AE()).data().begin();
	    index_type ld[N_];
	    for (XRANGE(i, N_)) { ld[i] = range::to()[i] - range::from()[i] - 1; }
	    operators_base::assign(p, ld);
	}

	size_type size1() const { return this->size(0); }
	size_type size2() const { return this->size(1); }
    private:
	const size_t N_;
    };

}

#endif /* _DISTRIBUTED_SLICE_HPP_ */
