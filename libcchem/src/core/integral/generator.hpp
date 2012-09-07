#ifndef INTEGRAL_GENERATOR_HPP
#define INTEGRAL_GENERATOR_HPP

#include "basis/basis.hpp"
#include "core/integral/aux.hpp"
#include "core/integral/evaluate.hpp"
#include <boost/tuple/tuple.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/intrinsic.hpp>
#include <boost/fusion/include/push_front.hpp>
#include <boost/fusion/include/vector_tie.hpp>

namespace integral {
namespace detail {

    template<size_t N, class E, class R>
    struct generator {
	typedef Basis::Shell Shell;
	typedef Basis::Block Block;
	// typedef const E& reference;
	// typedef const Shell& shell_ref;
	// typedef const Block& block_ref;

	typedef boost::array<R, N> array;
    
	const Basis &basis() const { return basis_; }
	const screening& screen() const { return screening_; }

	template<size_t I>
	const E& get() const { return array_.at(I); }
    
	const std::vector<Shell>& shells() const { return basis_.shells(); }

	generator<N+1,E,R> push_front(const E &e) const {
	    return generator<N+1,E,R>(*this, push_front(e, array_));
	}

	const Shell::Data& max_shell() const { return basis_.max(); }
	// const Basis::Block& max_block() const { return basis_.max_block(); }
	bool test(double value) const { return screening_.test(value); }
	bool test(const Shell &A, const Shell &B,
		  const Shell &C, const Shell &D) const {
	    return screening_.test(basis_.index(A), basis_.index(B),
				   basis_.index(C), basis_.index(D));
	}

	template<size_t N_, class E_, class R_>
	generator(const generator<N_,E_,R_> &other, const array &t)
	    : basis_(other.basis()), screening_(other.screen()), array_(t) {}
	generator(const Basis &basis, const screening &screening, const array &t)
	    : basis_(basis), screening_(screening), array_(t) {}

    private:
	const Basis &basis_;
	const screening &screening_;
	array array_;

	typename generator<2,E,R>::array
	push_front(const E &e, const typename generator<1,E,R>::array &A) const {
	    typename generator<2,E,R>::array array = {{ R(e), A.at(0) }};
	    return array;
	}

	typename generator<3,E,R>::array
	push_front(const E &e, const typename generator<2,E,R>::array &A) const {
	    typename generator<3,E,R>::array array = {{ R(e), A.at(0), A.at(1) }};
	    return array;
	}

	// typedef boost::reference_wrapper<const Shell> shell_ref;
	// static boost::array<shell_ref,2>
	// push_front(const Shell &A, const boost::array<shell_ref,1> &array) {
	// 	boost::array<shell_ref,2> a = {{ shell_ref(A), array[0] }};
	// 	return a;
	// }

    };
}

    template<size_t N>
    struct generator
	: detail::generator<N, Basis::Shell,
			    boost::reference_wrapper<const Basis::Shell> >
    {
	typedef Basis::Shell Shell;
	typedef boost::reference_wrapper<const Basis::Shell> R;
	typedef detail::generator<N, Shell, R> base_type;
	explicit generator(const base_type &base) : base_type(base) {}
	template<size_t N_, class E, class R>
	generator(const detail::generator<N_, E, R> &generator,
		  const typename base_type::array &data)
	    : base_type(generator, data) {}
	generator<N+1> push_front(const Shell &e) const {
	    return generator<N+1>(base_type::push_front(e));
	}
    };

    template<size_t N>
    struct block_generator
	: detail::generator<N, Basis::Block, Basis::Block>
    {
	typedef detail::generator<N, Basis::Block, Basis::Block> base_type;
	template<size_t N_>
	block_generator(const generator<N_> &generator,
			const typename base_type::array &data)
	    : base_type(generator, data) {}
	explicit block_generator(const base_type &base) : base_type(base) {}
    };

    // static generator<0>
    // make_generator(const Basis &basis, const screening &screening) {
    // 	generator<0>::array_type array;
    // 	return generator<0>(basis, screening, array);
    // }


    inline generator<1>
    make_generator(const Basis &basis, const screening &screening,
		   size_t i) {
	// using boost::fusion::make_vector;
	using boost::cref;
	typedef generator<1> G;
	G::array array = {{ cref(basis.shell(i)) }};
	return G(G::base_type(basis, screening, array));
    }

    inline generator<2>
    make_generator(const Basis &basis, const screening &screening,
		   size_t i, size_t j) {
	using boost::cref;
	typedef generator<2> G;
	G::array array = {{ cref(basis.shell(i)),
			    cref(basis.shell(j)) }};
	return G(G::base_type(basis, screening, array));
    }
		

}

#endif // INTEGRAL_GENERATOR_HPP
