#ifndef BOOST_NUMERIC_BINDINGS_CUBLAS_MATRIX_HPP
#define BOOST_NUMERIC_BINDINGS_CUBLAS_MATRIX_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <boost/numeric/bindings/cublas/forward.hpp>
#include <boost/numeric/bindings/cublas/assign.hpp>
#include <boost/numeric/bindings/cublas/vector.hpp>
#include <boost/numeric/bindings/cublas/storage.hpp>

#include <boost/numeric/bindings/begin.hpp>
#include <boost/numeric/bindings/end.hpp>
#include <boost/numeric/bindings/detail/adaptor.hpp>
#include <boost/numeric/bindings/detail/if_row_major.hpp>
#include <boost/numeric/bindings/detail/offset.hpp>
#include <boost/numeric/bindings/ublas/detail/convert_to.hpp>
#include <boost/numeric/bindings/ublas/storage.hpp>
#include <boost/numeric/bindings/ublas/matrix_expression.hpp>

#include <boost/assert.hpp>
#include <boost/static_assert.hpp>

namespace boost {
namespace numeric {
namespace bindings {
namespace cublas {

    struct matrix_tag;

    template<class E>
    struct matrix_expression : protected ublas::matrix_expression<E> {
	friend class ublas::matrix_expression<E>;
	// typedef cublas::matrix_expression<typename E::ublas_type> ublas_type;
	typedef ublas::matrix_expression<E> ublas_expression;
	typedef typename ublas_expression::expression_type expression_type;
	expression_type& operator()() {
	    return ublas_expression::operator()();
	}
	const expression_type& operator()() const {
	    return ublas_expression::operator()();
	}
    };


    template<class M, class U>
    struct matrix_base : protected U {
	//template< class, class> friend class traits::matrix_detail_traits;
	// template< class> friend class cublas::matrix_row;
	// template< class> friend class cublas::matrix_column;


	typedef U ublas_base_type;
	typedef cublas::matrix_tag matrix_tag;

	typedef U ublas_base;
	typedef typename U::size_type size_type;
	typedef typename U::difference_type difference_type;
	typedef typename U::value_type value_type;
	typedef typename U::reference reference;
	typedef typename U::const_reference const_reference;

	typedef typename U::storage_category storage_category;
	typedef typename U::orientation_category orientation_category;

	typedef typename U::iterator1 iterator1;
	typedef typename U::iterator2 iterator2;
	typedef typename U::const_iterator1 const_iterator1;
	typedef typename U::const_iterator2 const_iterator2;

	friend class bindings::detail::begin_impl<M, bindings::tag::value>;
	friend class bindings::detail::begin_impl<const M, bindings::tag::value>;
	template<class> friend class bindings::detail::get_dispatch;


	matrix_base() : U() {}
	template<class A0>
	matrix_base(A0 &a0) : U(a0) {}
	template<class A0, class A1>
	matrix_base(A0 &a0, A1 &a1) : U(a0, a1) {}
	template<class A0, class A1, class A2>
	matrix_base(A0 &a0, A1 &a1, const A2 &a2) : U(a0, a1, a2) {}

	using ublas_base::data;
	using ublas_base::size1;
	using ublas_base::size2;
	using ublas_base::find1;
	using ublas_base::find2;
	using ublas_base::begin1;
	using ublas_base::begin2;

	reference operator()(size_type i, size_type j) {
	    assert(false);
	    throw;
	}
	const_reference operator()(size_type i, size_type j) const {
	    assert(false);
	    throw;
	}

    };

    template<typename T, class A>
    struct matrix
	: cublas::matrix_expression<matrix<T,A> >,
	matrix_base<matrix<T,A>, ublas::matrix<T, ublas::column_major, A> >
    {
	typedef ublas::matrix<T, ublas::column_major, A> ublas_type;

	typedef matrix<T,A> self_type;
	typedef const self_type const_self_type;
	typedef ublas::matrix_reference<self_type> closure_type;
	typedef const ublas::matrix_reference<const_self_type> const_closure_type;

	typedef typename cublas::matrix_expression<matrix<T,A> >::expression_type expression_type;
	

	typedef matrix_base<self_type, ublas_type> base;

	typedef typename base::reference reference;
	typedef typename base::const_reference const_reference;
	typedef typename base::array_type array_type;
	typedef typename base::size_type size_type;

	typedef typename base::vector_temporary_type vector_temporary_type;
	typedef ublas::column_major functor_type;

	matrix() {
	    resize(0, 0, false);
	}
	matrix(const matrix& m) {
	    resize_assign(m);
	}
	matrix(size_type size1, size_type size2) {
	    resize(size1, size2, false);
	}
	matrix(size_type size1, size_type size2, const array_type &a)
	    : base(size1, size2, a) {}
	template<class E>
	matrix(const ublas::matrix_expression<E> &e) {
	    resize_assign(e());
	}
	template<class E>
	matrix(const cublas::matrix_expression<E> &e) {
	    resize_assign(e());
	}

	self_type operator=(const matrix& m) {
	    resize_assign(m);
	    return *this;
	}
	template<class E>
	self_type& operator=(const ublas::matrix_expression<E> &e) {
	    resize_assign(e());
	    return *this;
	}
	template<class E>
	self_type& operator=(const cublas::matrix_expression<E> &e) {
	    resize_assign(e);
	    return *this;
	}

	template<class E>
	matrix& operator+=(const cublas::matrix_expression<E> &e) {
	    matrix_plus_assign(*this, e);
	    return *this;
	}

	template<class M>
	void assign(const M &m) {
	    matrix_assign(*this, m);
	}

	void clear() { matrix_clear(*this); }

	void resize(size_type size1, size_type size2, bool preserve) {
	    BOOST_VERIFY(!preserve);
	    if (this->size1() == size1 && this->size2() == size2) return;
	    // std::cout << "resize: "
	    // 	      << this->size1() << " " << this->size2() << " "
	    // 	      << size1 << " " << size2 << std::endl;
	    matrix m(size1, size2, array_type(0));
	    this->swap(m);
	    this->data().swap(m.data());
	    size_t size = ublas::column_major::storage_size(size1, size2);
	    if (this->data().size() < size) {
		this->data().resize(size);
	    }
	}
    private:
	template<class M>
	void resize_assign(const M &m) {
	    resize(m.size1(), m.size2(), false);
	    assign(m);
	}
    };

    template<typename T, class L>
    struct matrix_adaptor
	: matrix_expression<matrix_adaptor<T,L> >,
	matrix_base<matrix_adaptor<T,L>, ublas::matrix<T, L, array_adaptor<T> > >
    {
	typedef matrix_adaptor<T,L> self_type;
	typedef const self_type const_self_type;
	typedef ublas::matrix_reference<self_type> closure_type;
	typedef const ublas::matrix_reference<const_self_type> const_closure_type;

	typedef array_adaptor<T> array_type;
	typedef typename array_adaptor<T>::pointer pointer;
	typedef matrix_base<
	    matrix_adaptor<T,L>, ublas::matrix<T, L, array_type> > base_type;

	// typedef matrix_adaptor closure_type;
	// typedef const matrix_adaptor const_closure_type;
	//typedef typename base_type::const_closure_type const_closure_type;

	matrix_adaptor(size_t size1, size_t size2, pointer data)
	    : base_type(size1, size2, array_type())
	{
	    size_t size = L::storage_size(size1, size2);
	    array_type array(size, data);
	    this->data().swap(array);
	} 
	void clear() {
	    matrix_clear(*this);
	}
    };

    template<class M, class U>
    struct matrix_vector
	: vector_expression<matrix_vector<M,U> >,
	  protected U
    {
	typedef U ublas_base;
	typedef typename ublas_base::value_type value_type;

	using ublas_base::size;
	using ublas_base::begin;
	using ublas_base::end;

	matrix_vector(M & data, size_t index)
	    : ublas_base(data, index) {}
	matrix_vector& operator=(const matrix_vector &e) {
	    assign(e);
	    return *this;
	}
	template<class E>
	matrix_vector& operator=(const vector_expression<E> &e) {
	    assign(e);
	    return *this;
	}
	template<class E>
	void assign(const vector_expression<E> &e) {
	    vector_assign(*this, e);
	}
    };

    template<class M>
    struct matrix_column: matrix_vector<M, ublas::matrix_row<M> > {
	typedef matrix_column self_type;
	typedef matrix_vector<M, ublas::matrix_row<M> > base;
	matrix_column(M & data, size_t index) : base(data, index) {}
	template<class E>
	self_type& operator=(const vector_expression<E> &e) {
	    base::operator=(e);
	    return *this;
	}
    };

    template<class M>
    struct matrix_row: matrix_vector<M, ublas::matrix_row<M> > {
	typedef matrix_row self_type;
	typedef matrix_vector<M, ublas::matrix_row<M> > base;
	matrix_row(M & data, size_t index) : base(data, index) {}
	template<class E>
	self_type& operator=(const vector_expression<E> &e) {
	    base::operator=(e);
	    return *this;
	}
    };



    template<class M>
    struct matrix_range
	: cublas::matrix_expression<matrix_range<M> >,
	matrix_base<matrix_range<M>, ublas::matrix_range<M> >
    {
	typedef ublas::matrix_range<M> ublas_base;
	typedef matrix_base<matrix_range<M>, ublas_base> base;

	matrix_range(const matrix_range& range)
	    : base(range) {}
	matrix_range(M &m, range r1, range r2)
	    : base(m, r1, r2) {}
	matrix_range& operator=(const matrix_range &m) {
	    matrix_assign(*this, m);
	    return *this;
	}
	template<class E>
	matrix_range& operator=(const cublas::matrix_expression<E>& m) {
	    matrix_assign(*this, m);
	    return *this;
	}
	template<class E>
	matrix_range& operator+=(const cublas::matrix_expression<E>& m) {
	    matrix_plus_assign(*this, m);
	    return *this;
	}
	template<class E>
	matrix_range& operator=(const ublas::matrix_expression<E> &e) {
	    assign(e());
	    return *this;
	}

	template<class E>
	void assign(const E &m) {
	    matrix_assign(*this, m);
	}
    };

    template<class M>
    matrix_range<M> project(M &m, range r1, range r2) {
	return matrix_range<M>(m, r1, r2);
    }

    template<class M>
    matrix_range<M> subrange(M &m,
			     size_t start1, size_t stop1,
			     size_t start2, size_t stop2) {
	return matrix_range<M>(m, range(start1, stop1), range(start2, stop2));
    }

}
}
}
}

namespace boost {
namespace numeric {
namespace bindings {
namespace detail {

template< class M, typename Id >
struct adaptor< M, Id, typename boost::enable_if<
			   boost::is_same<
			       typename M::matrix_tag,
			       cublas::matrix_tag>
			   >::type>
: adaptor<typename M::ublas_base_type,
	  typename copy_const<Id, typename M::ublas_base_type>::type >
{};

}
}
}
}

#endif // BOOST_NUMERIC_BINDINGS_CUBLAS_MATRIX_HPP
