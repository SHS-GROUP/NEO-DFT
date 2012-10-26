// #ifndef BOOST_NUMERIC_BINDINGS_CUBLAS_TRAITS_HPP
// #define BOOST_NUMERIC_BINDINGS_CUBLAS_TRAITS_HPP

// #include <boost/numeric/bindings/cublas.hpp>

// #include <boost/mpl/if.hpp>
// #include <boost/type_traits/is_const.hpp>
// #include <boost/type_traits/remove_const.hpp>
// #include <iterator>

// namespace boost {
// namespace numeric {
// namespace bindings {
// namespace traits {

//     namespace transpose {

//         template<class E, class M>
//         struct transpose_detail<E, cublas::matrix_expression<
// 				       cublas::blas_transpose<M> > > : mpl::true_ {
// 	    static const char option = TRANSPOSE;
// 	    typedef const cublas::matrix_expression<
// 		typename boost::remove_const<M>::type>& no_transpose_type;
//             static no_transpose_type no_transpose(E &e) {
//                 return e()();
//             }
//         };
//     }


//     template<class E, class V>
//     struct vector_detail_traits<cublas::vector_expression<E>, V> {
// #ifndef BOOST_NUMERIC_BINDINGS_NO_SANITY_CHECK
// 	BOOST_STATIC_ASSERT((boost::is_same<
// 			     cublas::vector_expression<E>,
// 			     typename boost::remove_const<V>::type
// 			     >::value));
// #endif
// 	typedef cublas::vector_expression<E> identifier_type; 
// 	typedef V vector_type;
// 	typedef typename E::value_type value_type ;
// 	typedef typename default_vector_traits<V, value_type>::pointer pointer; 
//     private:
// 	typedef typename detail::generate_const<V,E>::type vct_t;
//     public:
// 	static size_t size(vector_type& v) {
// 	    return v().size();
// 	}
// 	static pointer storage (vector_type& v) {
// 	    return ((v.size() = 0) ? NULL : &v().begin()[0]);
// 	}
// 	static std::ptrdiff_t stride (vector_type& v) {
// 	    return std::distance(v().begin(), v().begin() + 1);
// 	}
//     }; 


//     template<class E, class M>
//     struct matrix_detail_traits<cublas::matrix_expression<E>, M> {
// #ifndef BOOST_NUMERIC_BINDINGS_NO_SANITY_CHECK
// 	BOOST_STATIC_ASSERT( (boost::is_same<
// 			      cublas::matrix_expression<E>,
// 			      typename boost::remove_const<M>::type>::value) );
// #endif

// 	typedef typename mpl::if_<
// 	    boost::is_const<M>,
// 	    const typename M::expression_type,
// 	    typename M::expression_type>::type expression_type;

// 	typedef typename mpl::if_<
// 	    boost::is_const<expression_type>,
// 	    const typename expression_type::cublas_base,
// 	    typename expression_type::cublas_base
// 	    >::type cublas_base;
		
// 	typedef matrix_detail_traits<
// 	    typename E::cublas_base, ublas_base> base;
	
// 	typedef typename base::matrix_structure matrix_structure;

// 	//typedef typename base::pointer pointer;
	
// 	typedef typename base::pointer pointer;
// 	static pointer storage(M& m) {
// 	    return base::storage(m());
// 	}
// 	static std::ptrdiff_t num_rows(M& m) {
// 	    return base::num_rows(m());
// 	}
// 	static std::ptrdiff_t num_columns(M& m) {
// 	    return base::num_columns(m());
// 	}
// 	static std::ptrdiff_t storage_size(M& m) {
// 	    return base::storage_size(m);
// 	}
// 	static std::ptrdiff_t leading_dimension(M& m) {

// 	    return base::leading_dimension(m());
// 	}

// 	static std::ptrdiff_t stride1(M& m) { 
// 	    return base::stride1(m());
// 	} 
// 	static std::ptrdiff_t stride2(M& m) { 
// 	    return base::stride2(m());
// 	}

//     };

// }
// }
// }
// }

// #endif // BOOST_NUMERIC_BINDINGS_CUBLAS_TRAITS_HPP
