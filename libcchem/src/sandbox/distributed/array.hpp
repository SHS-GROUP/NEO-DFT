#ifndef _DISTRIBUTED_ARRAY_HPP_
#define _DISTRIBUTED_ARRAY_HPP_


#include <ga++.h>
#include <boost/mpl/if.hpp>
#include <boost/type_traits.hpp>
#include <boost/static_assert.hpp>

#define static_assert BOOST_STATIC_ASSERT

#include "util/array.hpp"
#include "distributed/slice.hpp"
#include "distributed/element.hpp"
#include "distributed/ops.hpp"
#include "distributed/traits.hpp"

namespace mpl = boost::mpl;


// 	template<class A>
// 	typename A::range_n operator()(A &a) {
// 	    static_assert(array_traits<A>::N == N);
// 	    return typename A::slice_type();
// 	}

namespace GA {

    template<typename T>
    struct type { static const int value = 0; };

    template<>
    struct type<double> { static const int value = MT_F_DBL; };

}

namespace distributed {


    template<typename T>
    struct scalar_operators {
	scalar_operators(GA::GlobalArray &ga) : ga(ga) {}
	void operator=(T t) { (t == T(0)) ? ga.zero() : ga.fill(&t); }
	void operator+=(T t) { ga.addConstant(&t); }
	void operator-=(T t) { ga.addConstant(&T(-t)); }
	void operator*=(T t) { ga.scale(&t); }
	void operator/=(T t) { ga.scale(&(T(1)/t)); }
    private:
	GA::GlobalArray &ga;
    };

    template<typename T>
    struct array_base : GA::GlobalArray, scalar_operators<T> {
	typedef GA::GlobalArray ga_base;
	typedef scalar_operators<T> scalar_base;
	typedef T value_type;
	typedef T* pointer;
	typedef const T* const_pointer;

	typedef int index_type;
	typedef size_t size_type;
	static const int ga_type = GA::type<T>::value;

    protected:
	template<size_t N>
	array_base(const index_type (&dims)[N])
	    : ga_base(ga_type, N, const_cast<index_type*>(dims), (char*)"", NULL),
	      scalar_base(static_cast<ga_base&>(*this)) { }

	void assign(T t) { scalar_base::operator=(t); }

	template<class E>
	void assign(const E &e) {
	    distributed::assign(*this, e);
	}
    };

    template<typename T, size_t N>
    struct array : array_base<T> {
	typedef array_base<T> base;
	typedef typename base::index_type index_type;
	typedef util::Array<index_type,N> index_array;
	typedef distributed::element<T,N> element;
	typedef distributed::const_element<T,N> const_element;

	array(index_array dims) : base(dims.elems) {}

	template<class E>
	void operator=(const E &e) {
	    typedef typename  mpl::if_<boost::is_arithmetic<E>, T, E>::type type;
	    this->assign(*this, type(e));
	}

	index_array dims() const  {
	    int type, ndim;
	    index_array dims;
	    this->inquire(&type, &ndim, dims);
	    return dims;
	}

 	template<size_t M>
	array_slice<T,M> operator[](const util::range_n<M> range) {
	    static_assert(M <= N);
	    typedef array_slice<T,M> S;
	    const index_array from = 0, to = this->dims() - 1;
	    return S(*this, range.lower(from), range.upper(to));
	}




    };

}

#endif /* _DISTRIBUTED_ARRAY_HPP_ */
