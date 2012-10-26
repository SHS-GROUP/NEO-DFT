#ifndef BOOST_PP_IS_ITERATING

#ifndef TENSOR_TENSOR_HPP
#define TENSOR_TENSOR_HPP

#include "tensor/config.hpp"
#include "tensor/forward.hpp"
#include "tensor/storage/storage.hpp"
#include "tensor/index.hpp"
#include "tensor/view.hpp"
#include "tensor/functional.hpp"
#include "tensor/operators.hpp"

#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/iteration/iterate.hpp>
#include <boost/preprocessor/enum_params.hpp>

#include <boost/fusion/include/transformation.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/mpl.hpp>
#include <boost/mpl/count_if.hpp>
#include <boost/type_traits.hpp>


namespace tensor {
namespace detail {


    template<class I>
    struct index_ranges {

	template<class T>
	struct remove_cr {
	    typedef typename boost::remove_const<
		typename boost::remove_reference<T>::type
		>::type type;
	};

	template<class T>
	struct is_index_ :
	    detail::is_index<typename remove_cr<T>::type> {};

	typedef typename boost::fusion::result_of::as_map<
	    typename boost::fusion::result_of::filter_if<
		typename boost::mpl::transform<
		    I, remove_cr<boost::mpl::_> >::type,
		is_index_<boost::mpl::_> >::type
	    >::type type;
	typedef typename boost::fusion::result_of::keys<type>::type keys_type;
	static const size_t rank = boost::mpl::size<type>::value;

	template<class T>
	struct transform {
	    template<class R>
	    void operator()(R&) const {}
	    template<int K, class R>
	    void operator()(index<K,R> &e) const {
		BOOST_AUTO(begin, e.begin());
		BOOST_AUTO(end, e.end());
		if (begin == e.min()) begin = 0;
		if (end == e.max()) end = data_.size(i_++);
		e = index<K,R>(begin, end, e.increment());
	    }
	    explicit transform(const T &data) : data_(data), i_(0) {}
	    const T &data_;
	    mutable int i_;
	};

	template<class T, class S>
	static type generate(const T &t, const S& s) {
 	    BOOST_AUTO(s_, fusion::as_vector(s));
	    fusion::for_each(s_, transform<T>(t));
	    return fusion::as_map(fusion::filter_if<is_index_<boost::mpl::_> >
	     			  (fusion::as_vector(s_)));
	}
    };

    template<class A, size_t N = detail::traits<A>::rank>
    struct reference_generator {
	typedef reference_generator self_type;
	//typedef detail::generator<self_type> generator;
	static const size_t rank = N;

	reference_generator(A &data) : data_(data) {}
    private:
	A &data_;

    public:
	template<class I>
	struct contains_index :
	    boost::mpl::count_if<
	    I, detail::is_index<
		   boost::remove_reference<boost::mpl::_> > >::type {};

	template<class F>
	struct result {
	    BOOST_MPL_ASSERT_MSG((false), MUST_NOT_INSTANTIATE, (F));
	};

	template<class F, class I>
	struct result<F(I)> {
	    typedef typename index_ranges<I>::type I_;
	    static const size_t rank = index_ranges<I>::rank;

	    template<class T, class U>
	    struct if_const {
		typedef typename boost::mpl::if_<
		    boost::is_const<F>, T, U>::type type;
	    };

	    typedef typename if_const<
		typename storage::const_array_view<A, rank>::type,
		typename storage::array_view<A, rank>::type
		>::type A_;

	    typedef typename boost::mpl::eval_if_c<
		rank,
		if_const<
		    const_tensor_view<A_, I_>,
		    tensor_view<A_, I_>
		    >,
		if_const<
		    typename detail::traits<A>::const_reference,
		    typename detail::traits<A>::reference
		    >
		>::type type;
	};

	// template<class G, class I>
	// struct result<self_type(G&, I)> : result<A(I)> {};

	// template<class G, class I>
	// struct result<const self_type(const G&, I)> : result<const A(I)> {};

	template<class G, class I>
	typename result<const self_type(I)>::type
	operator()(const G&, const I &indices) const {
	    typedef result<const self_type(I)> R;
	    typedef storage::const_array_view<A, R::rank> V;
	    return typename R::type(V::generate(data_, indices));
	}

	template<class G, class I>
	typename result<self_type(I)>::type
	operator()(const G&, const I &indices) {
	    typedef result<self_type(I)> R;
	    typedef storage::array_view<A, R::rank> V;
	    return typename R::type(V::generate(data_, indices));
	}

    // private:
    // 	template<class T, size_t Rank, class I>
    // 	static T generate(boost::mpl::size_t<Rank>, A &data, const I &indices) {
    // 	    return T(storage::array_view<A,Rank>::generate(data, indices));
    // 	}
    // 	template<class T, size_t Rank, class I>
    // 	static T generate(boost::mpl::size_t<Rank>, const A &data, const I &indices) {
    // 	    return T(storage::const_array_view<A,Rank>::generate(data, indices));
    // 	}
    };


    template<class A>
    struct const_tensor_base
	: const_generator<reference_generator<A> >
    {
    private:
	typedef const_generator<reference_generator<A> > generator;
	typedef typename detail::traits<A> traits;

    public:
	typedef const_tensor_base self_type;
	static const size_t rank = traits::rank;
	// typedef typename boost::mpl::print<same_rank<A,N> >::type print;
	typedef typename traits::value_type value_type;
	typedef typename traits::const_reference const_result_type;

	const_tensor_base(const A &data) :
	    generator(reference_generator<A>(data_)),
	    data_(data) {}

	const A& data() const { return data_; }

	template<typename I>
	const_result_type operator[](const I &indices) const;
	
	template<class O>
	struct result : generator::template result<O> {};

	using generator::operator();

	size_t size(size_t i) const {
	    return storage::size(data_, i);
	}
    protected:
	A data_;
    };


    template<class A>
    struct tensor_base :
	const_tensor_base<A>,
	detail::generator<reference_generator<A> >
    {
    private:
	typedef tensor_base self_type;
	typedef const_tensor_base<A> base_type;
	typedef detail::generator<reference_generator<A> > generator;

    public:
	static const size_t rank = base_type::rank;
	typedef typename detail::traits<A>::reference result_type;
	typedef typename detail::traits<A>::reference reference;

	tensor_base(const A &data = A()) :
	    base_type(data),
	    generator(reference_generator<A>(base_type::data_)) {}

	A& data() { return this->data_; }


	template<class O>
	struct result : generator::template result<O> {};

	using base_type::operator[];

	template<typename I>
	result_type operator[](const I &indices);
	
	typename boost::mpl::if_c<
	    (rank == 1),
	    reference,
	    tensor_ref<typename storage::array_ref<A>::type>
	    >::type			  
	operator[](const int &i) {
	    typedef typename boost::mpl::if_c<
		(rank == 1),
		reference,
		tensor_ref<typename storage::array_ref<A>::type>
		    >::type T;
	    return T(storage::array_ref<A>::generate(this->data(), i));
	}

	using base_type::operator();
	using generator::operator();

    };


}
}


namespace tensor {

    template<class A>
    struct tensor_ref :
	detail::tensor_base<A>,
	expression<tensor_ref<A> >
    {
	typedef detail::tensor_base<A> base_type;

	using base_type::operator();
	using base_type::operator[];

	explicit tensor_ref(const A &data) : base_type(data) {}
    };

}

// generate specializations
#define BOOST_PP_ITERATION_LIMITS (1, TENSOR_MAX_RANK)
#define BOOST_PP_FILENAME_1       "tensor/tensor.hpp" // this file
#include BOOST_PP_ITERATE()

#endif // TENSOR_TENSOR_HPP

#else // BOOST_PP_IS_ITERATING

#define N BOOST_PP_ITERATION()


namespace tensor {

    template<typename T, class A>
    struct Tensor<N, T, A> :
	detail::tensor_base<A>,
	expression<Tensor<N, T, A> >
    {
    	typedef detail::tensor_base<A> base_type;
	typedef size_t size_type;

	using base_type::operator();
	using base_type::operator[];

	Tensor() {}
	explicit Tensor(BOOST_PP_ENUM_PARAMS(N, size_type size)) {
	    size_type dims[] = { BOOST_PP_ENUM_PARAMS(N, size) };
	    resize(dims);
	}
	explicit Tensor(const size_type (&dims)[N]) {
	    resize(dims);
	}
	void resize(const size_type (&dims)[N]) {
	    storage::array<A>::resize(this->data(), dims);
	}
    };

}

#undef N

#endif // BOOST_PP_IS_ITERATING
