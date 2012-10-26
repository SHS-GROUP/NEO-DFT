#ifndef BOOST_PP_IS_ITERATING

#ifndef TENSOR_GENERATOR_HPP
#define TENSOR_GENERATOR_HPP

#include "tensor/config.hpp"
#include "tensor/index.hpp"

#include <boost/mpl/assert.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/vector_tie.hpp>

#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/iteration/iterate.hpp>
#include <boost/preprocessor/enum_params.hpp>


namespace tensor {
namespace detail {

    template<class G, size_t N = G::rank>
    struct generator;

    template<class G, size_t N = G::rank>
    struct const_generator;        


    template<class T>
    struct fail_instantiate {
	BOOST_MPL_ASSERT_MSG((false), FAIL_INSTANTIATE, (T));
    };

}
}


// generate specializations
#define BOOST_PP_ITERATION_LIMITS (1, TENSOR_MAX_RANK)
#define BOOST_PP_FILENAME_1       "tensor/generator.hpp" // this file
#include BOOST_PP_ITERATE()

#endif // TENSOR_GENERATOR_HPP

#else // BOOST_PP_IS_ITERATING

#define N BOOST_PP_ITERATION()


namespace tensor {
namespace detail {

    template<class G>
    struct const_generator<G,N> {
	const_generator(const G &g = G()) : g_(g) {}

	template<class O>
	struct result : G::template result<O> {};

	template<BOOST_PP_ENUM_PARAMS(N, class I)>
	typename result<
	    const G(indices<boost::fusion::vector<BOOST_PP_ENUM_PARAMS(N, I)> >)
	    >::type
	operator()(BOOST_PP_ENUM_BINARY_PARAMS(N, const I, &i)) const {
	    using boost::fusion::vector_tie;
	    return g_(*this, vector_tie(BOOST_PP_ENUM_PARAMS(N, i)));
	}

	template<class I>
	typename result<const G(indices<I>)>::type
	operator()(const indices<I> &indices) const {
	    return g_(*this, indices);
	}

    private:
	G g_;
    };

    template<class G>
    struct generator<G,N> {
	generator(const G &g = G()) : g_(g) {}

	template<class O>
	struct result : G::template result<O> {};

	template<BOOST_PP_ENUM_PARAMS(N, class I)>
	typename result<
	    G(indices<boost::fusion::vector<BOOST_PP_ENUM_PARAMS(N, I)> >)
	    >::type
	operator()(BOOST_PP_ENUM_BINARY_PARAMS(N, const I, &i)) {
	    using boost::fusion::vector_tie;
	    return g_(*this, vector_tie(BOOST_PP_ENUM_PARAMS(N, i)));
	}

	template<class I>
	typename result<G(indices<I>)>::type
	operator()(const indices<I> &indices) {
	    return g_(*this, indices);
	}

	// template<class I>
	// typename result<const G(indices<I>)>::type
	// operator()(const indices<I> &indices) const {
	//     return g_(*this, indices);
	// }

    private:
	G g_;
    };


} // namespace detail
} // namespace tensor

#undef N

#endif // BOOST_PP_IS_ITERATING
