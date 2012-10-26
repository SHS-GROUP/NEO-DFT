#ifndef TENSOR_STORAGE_MULTI_ARRAY_HPP
#define TENSOR_STORAGE_MULTI_ARRAY_HPP

#include "tensor/index.hpp"
#include "tensor/lambda.hpp"

#include <boost/multi_array.hpp>

#include <boost/mpl/identity.hpp>
#include <boost/mpl/contains.hpp>
#include <boost/mpl/count_if.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/assert.hpp>

#include <boost/fusion/include/transformation.hpp>
#include <boost/fusion/include/intrinsic.hpp>
#include <boost/fusion/include/map.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/mpl.hpp>
#include <boost/fusion/include/boost_array.hpp>
#include <boost/fusion/include/io.hpp>
#include <boost/type_traits.hpp>


namespace tensor {
namespace storage {

    //namespace mpl = boost::mpl;

    template<class A, size_t N = A::dimensionality>
    struct multi_array_indices {

	typedef typename A::index index;
	typedef typename A::size_type size_type;
	typedef typename A::index_range index_range;

	template<class R>
	struct tuple {
	    static const size_t size = boost::mpl::count_if<
		R, detail::is_index<
		       boost::remove_const<
			   boost::remove_reference<boost::mpl::_> > >
		>::value;
	typedef typename boost::mpl::if_c<
	    (size > 0),
		typename A::index_gen::template gen_type<size, N>::type,
		boost::array<index, N>
		>::type type;
	};

	template<class V, class R>
	static
	typename boost::enable_if_c<(tuple<R>::size == 0), V>::type
	generate_view(A &a, const R &r) {
	    typename tuple<R>::type indices;
	    BOOST_AUTO(it, indices.rbegin());
	    boost::fusion::for_each(r, (*lambda::ref(it)++ = lambda::_1));
	    return a(indices);
	}

	template<class V, class R>
	static
	typename boost::enable_if_c<(tuple<R>::size > 0), V>::type
	generate_view(A &a, const R &r) {
	    typedef typename tuple<R>::type I;
	    return a[generate<I>(a.index_bases(), a.shape(),
				 r, boost::indices)];
	}

	template<class T, class R, class G>
	static typename boost::disable_if<boost::mpl::empty<R>, T>::type
	generate(const index *bases, const size_type *shape,
		 const R &r, const G &g) {
	    namespace fusion = boost::fusion;
	    BOOST_AUTO(const &gr,
		       g[as_range(*bases++, *shape++, fusion::back(r))]);
	    return generate<T>(bases, shape, fusion::pop_back(r), gr);
	}

	template<class T, class R, class G>
	static typename boost::enable_if<boost::mpl::empty<R>, T>::type
	generate(const index*, const size_type*,
		 const R &r, const G &g) {
	    return T(g);
	}
	    
	template<class R>
	static const R& 
	as_range(const index &base, const size_type &size, const R &i) {
	    return i;
	}

	template<int K, class R>
	static index_range
	as_range(const index &base, const size_type &size,
		 const ::tensor::index<K,R> &i) {
	    typedef typename index_range::index index;
	    typedef typename index_range::size_type size_type;
	    return (i.all() ?
		    index_range(base, base+size) :
		    index_range(*i.begin(), *i.end(), i.increment()));
	}

    };


    template<class A>
    struct multi_array_impl_base {
	typedef void enable;

	typedef A type;
	typedef typename A::size_type size_type;
	typedef typename A::index stride_type;
	typedef typename A::index index;
	typedef typename A::index_range index_range;

	typedef typename A::reference array_ref;
	typedef typename A::const_reference const_array_ref;

	static const size_t rank = A::dimensionality;
	typedef typename A::template array_view<1>::type S;
	typedef typename S::value_type value_type;
	typedef typename S::reference reference;
	typedef typename S::const_reference const_reference;

	typedef value_type* pointer;
	typedef const value_type* const_pointer;

	typedef typename S::iterator iterator;
	typedef typename S::const_iterator const_iterator;


	template<size_t N>
	struct const_view {
	    typedef typename boost::mpl::eval_if_c<
		(N > 0),
		    typename A::template const_array_view<N>,
		    boost::mpl::identity<const_reference>
		    >::type type;
	};

	template<size_t N>
	struct view {
	    typedef typename boost::mpl::eval_if_c<
		(N > 0),
		    typename A::template array_view<N>,
		    boost::mpl::identity<reference>
		    >::type type;
	};
	
	template<class V, class A_, class R>
	static V generate_view(A_ &a, const R &r) {
	    return multi_array_indices<A_>::template generate_view<V>(a, r);
	}
	
	template<class K>
	struct push_back_if_in {
	    template<class F> struct result;
	    
	    template<class F, class S, class E>
	    struct result<F(S&, E&)> {
		typedef typename boost::fusion::result_of::first<E>::type K_;
		typedef typename boost::fusion::result_of::second<E>::type E_;
		typedef typename boost::mpl::eval_if<
		    boost::mpl::contains<K,K_>,
		    boost::fusion::result_of::push_back<S,E_>,
		    boost::add_reference<S> >::type type;

	    static type apply(S& s, E& e) {
		return apply(s, e, boost::mpl::contains<K,K_>());
	    }

	    template<class S_>
	    static type apply(const S_& s, E& e, boost::mpl::false_) {
		return s;
	    }

	    template<class S_>
	    static type apply(const S_& s, E& e, boost::mpl::true_) {
		return boost::fusion::push_back<S_,E_>(s, e.second);
	    }

	};
	    
	template<class S, class E>
	typename result<push_back_if_in(const S&, const E&)>::type
	operator()(const S& s, const E& e) const {
	    typedef result<push_back_if_in(const S&, const E&)> R;
	    return R::apply(s, e);
	}

	};


	template<typename R, class I>
	static R generate(A &a, const I &i) {
	    return generate(boost::type<R>(), a, i);
	}

	template<class I>
	static const_array_ref
	generate(boost::type<const_array_ref>, A &a, const I &i) {
	    return a[i];
	}

	template<class I>
	static array_ref
	generate(boost::type<array_ref>, A &a, const I &i) {
	    return a[i];
	}

	template<typename R, class A_, class K, class M>
	static R element(A_ &a, const K& keys, const M& map) {
	    namespace fusion = boost::fusion;
	    BOOST_AUTO(const &seq,
		       fusion::fold(fusion::as_map(map),
				    fusion::vector<>(),
				    push_back_if_in<K>()));

	    typedef BOOST_TYPEOF(seq) S;
	    static const size_t N  = boost::mpl::size<S>::value;
	    BOOST_MPL_ASSERT_MSG((N == rank),
				 TENSOR_INDEX_INCORRECT_KEYS,
				 (K, M, S));
				  
	    boost::array<typename A_::index, rank> indices;
	    BOOST_AUTO(it, indices.rbegin());
	    fusion::for_each(seq, (*lambda::ref(it)++ = lambda::_1));
	    return a(indices);
	}

	static pointer data(A &a) {
	    boost::array<index, rank> indices;
	    std::copy(a.index_bases(), a.index_bases()+rank, indices.begin());
	    return &a(indices);
	}

	static iterator begin(A &a) {
	    return a.begin();
	}

	static iterator end(A &a) {
	    return a.end();
	}

	static size_type size(const A &a, int dimension) {
	    return a.shape()[rank-(dimension+1)];
	}

	static stride_type stride(const A &a, int dimension) {
	    return a.strides()[rank-(dimension+1)];
	}

    };


    template<class A>
    struct multi_array_impl {};

    template<class A>
    struct multi_array_impl<const A>
	: multi_array_impl<A> {};

    template<typename T, size_t N, class A>
    struct multi_array_impl<boost::multi_array<T,N,A> >
	: multi_array_impl_base<boost::multi_array<T,N,A> >
    {
	typedef boost::multi_array<T,N,A> array_type;
	template<typename S>
	static void resize(boost::multi_array<T,N,A> &a, const S (&dims)[N]) {
	    boost::array<S,N> dims_;
	    std::copy(dims, dims + N, dims_.rbegin());
	    a.resize(dims_);
	}
    };

    template<typename T, size_t N>
    struct multi_array_impl<
	boost::detail::multi_array::multi_array_view<T, N> > :
    multi_array_impl_base<
	boost::detail::multi_array::multi_array_view<T, N> > {};

    template<typename T, size_t N, class U>
    struct multi_array_impl<
	boost::detail::multi_array::const_multi_array_view<T, N, U*> > :
    multi_array_impl_base<
	boost::detail::multi_array::const_multi_array_view<T, N, U*> > {};

    template<typename T, size_t N>
    struct multi_array_impl<
	boost::detail::multi_array::sub_array<T, N> > :
    multi_array_impl_base<
	boost::detail::multi_array::sub_array<T, N> > {};

    template<typename T, size_t N, class U>
    struct multi_array_impl<
	boost::detail::multi_array::const_sub_array<T, N, U*> > :
    multi_array_impl_base<
	boost::detail::multi_array::const_sub_array<T, N, U*> > {};

    template<class A>
    struct array_traits<A, typename multi_array_impl<A>::enable>
	: multi_array_impl<A> {};



} // namespace storage
} // namespace tensor

#endif // TENSOR_STORAGE_MULTI_ARRAY_HPP
