#ifndef BOOST_BINDINGS_CUBLAS_EXPRESSION_HPP
#define BOOST_BINDINGS_CUBLAS_EXPRESSION_HPP

#include <boost/numeric/bindings/cublas/forward.hpp>
//#include <boost/numeric/bindings/cublas/cublas3.hpp>

// #include <boost/numeric/bindings/begin.hpp>
// #include <boost/numeric/bindings/end.hpp>
#include <boost/numeric/bindings/detail/adaptor.hpp>
// #include <boost/numeric/bindings/detail/if_row_major.hpp>
// #include <boost/numeric/bindings/detail/offset.hpp>
// #include <boost/numeric/bindings/ublas/detail/convert_to.hpp>
// #include <boost/numeric/bindings/ublas/storage.hpp>
// #include <boost/numeric/bindings/ublas/matrix_expression.hpp>

namespace boost {
namespace numeric {
namespace bindings {
namespace cublas {

    template<class E>
    struct transpose : matrix_expression<transpose<E> > {
	//typedef typename E::const_closure_type expression_closure_type;
	explicit transpose(E &e) : e_(e) {}
	const E& expression() const { return e_; } 
	size_t size1() const { return e_.size1(); }
	size_t size2() const { return e_.size2(); }
    private:
	E &e_;
    };

    template<class E>
    struct matrix_expression<transpose<E> > {
	const transpose<E>& operator()() const {
	    return *static_cast<const cublas_transpose<E>*>(this);
	}
    };

    template<class E>
    transpose<E> trans(E &e) {
	return transpose<E>(e);
    }
    

    template<class E>
    struct cublas_matrix_expression
	: matrix_expression<cublas_matrix_expression<E> >
    {
	typedef E expression_type;
	const expression_type& operator()() const {
	    return *static_cast<const expression_type*>(this);
	}
	size_t size1() const { return operator()().size1(); }
	size_t size2() const { return operator()().size2(); }
	template<typename T, class C>
	void operator()(T beta, cublas::matrix_expression<C> &c) const {
	    operator()()(beta, c);
	}
    };

    template<class E>
    struct matrix_expression<cublas_matrix_expression<E> > {
	typedef cublas_matrix_expression<E> expression_type;
	const expression_type& operator()() const {
	    return *static_cast<const expression_type*>(this);
	}
	size_t size1() const { return operator()().size1(); }
	size_t size2() const { return operator()().size2(); }
    };


    template<class A, class B>
    void matrix_assign(cublas::matrix_expression<A> &a,
		       const cublas::matrix_expression<
		       cublas_matrix_expression<B> > &b) {
	b()(typename A::value_type(0), a);
    }

    template<class A, class B>
    void matrix_plus_assign(cublas::matrix_expression<A> &a,
			    const cublas::matrix_expression<
			    cublas_matrix_expression<B> > &b) {
	b()(typename A::value_type(1), a);
    }


    template<class A, class B>
    struct gemm_expression : cublas_matrix_expression<gemm_expression<A,B> > {
	gemm_expression(const A &a, const B &b) : a_(a), b_(b) {}
	size_t size1() const { return a_.size1(); }
	size_t size2() const { return b_.size2(); }
	template<typename T, class C>
	void operator()(T beta, cublas::matrix_expression<C> &c) const {
	    gemm(T(1), a_, b_, beta, c);
	}
    private:
	const A &a_;
	const B &b_;
    };

    template<class A, class B>
    gemm_expression<A, B>
    gemm(const A &a, const B &b) {
      // std::cout << a.size1() << std::endl;
	return gemm_expression<A, B>(a, b);
    }

}
}
}
}


namespace boost {
namespace numeric {
namespace bindings {
namespace detail {

template< class M, typename Id, class Enable >
struct adaptor< cublas::transpose<M>, Id, Enable >
{

    typedef typename copy_const< Id, M >::type adapted_type;
    typedef typename property_map_of< adapted_type >::type map;

    typedef mpl::map<
        mpl::pair<tag::value_type, typename mpl::at<map, tag::value_type>::type>,
        mpl::pair<tag::entity, typename mpl::at<map, tag::entity>::type>,
        mpl::pair<tag::size_type<1>, typename mpl::at<map, tag::size_type<1> >::type>,
	mpl::pair<tag::size_type<2>, typename mpl::at<map, tag::size_type<2> >::type>,
	mpl::pair<tag::data_structure, typename mpl::at<map, tag::data_structure>::type>,

        mpl::pair<tag::data_order,
		  typename mpl::if_<
		      is_same<
			  typename mpl::at<map, tag::data_order>::type,
			  tag::row_major>,
		      tag::column_major,
		      tag::row_major
		      >::type>,

        mpl::pair<tag::stride_type<1>,
		  typename mpl::at<map, tag::stride_type<2> >:: type>,
        mpl::pair<tag::stride_type<2>,
		  typename mpl::at<map, tag::stride_type<1> >:: type>
	> property_map;

    static std::ptrdiff_t size1( const Id& id ) {
        return id.size1();
    }

    static std::ptrdiff_t size2( const Id& id ) {
        return id.size2();
    }

    static typename result_of::begin_value< adapted_type >::type
    begin_value( Id& id ) {
        return bindings::begin_value( id.expression() );
    }

    static typename result_of::end_value< adapted_type >::type
    end_value( Id& id ) {
        return bindings::end_value( id.expression() );
    }

    static std::ptrdiff_t stride1( const Id& id ) {
        return bindings::stride2( id.expression() );
    }

    static std::ptrdiff_t stride2( const Id& id ) {
        return bindings::stride1( id.expression() );
    }
};

}
}
}
}


#endif // BOOST_BINDINGS_CUBLAS_EXPRESSION_HPP
