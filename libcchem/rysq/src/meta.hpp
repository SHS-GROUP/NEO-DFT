#ifndef RYSQ_META_HPP
#define RYSQ_META_HPP

#include <boost/mpl/min_max.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/int.hpp>

#include "rysq-core.hpp"
#include "transpose.hpp"

namespace rysq {
    namespace meta {

	namespace mpl = boost::mpl;

	template<int n> struct binomial2 {
	    static const int value = (n*n - n)/2;
	};

	template<int x, int y> struct max {
	    static const int value = (x >= y) ? x : y;
	};

	template<int x, int y> struct min {
	    static const int value = (x < y) ? x : y;
	};

	template<rysq::type S>
	struct L : boost::mpl::int_<(S < 0) ? -S : S> {};

	template<rysq::type S> class shell {
	public:
	    static const rysq::type type = S;
	    static const int L = (S < 0) ? -S : S;

	    static const int nc = (S < 0) ? -S + 1 : 1;
	    static const int size = binomial2<L+2>::value + (nc-1);
	    static const int index = (S < 0) ? 0 : (S*(S + 1)*(S + 2))/6;

	    typedef mpl::integral_c_tag tag;
	    static const int value = ((1 << L) | ((S < 0) ? ((1 << L) - 1) : 0)) - 1;
	};

	template<rysq::type A_, rysq::type B_>
	class state {
	public:
	    typedef shell<A_> A;
	    typedef shell<B_> B;
	    static const int L = A::L + B::L;
	    static const int nc0 = A::nc;
	    static const int nc1 = B::nc;
	    static const int nc = nc0*nc1;
	    static const int size = A::size*B::size;

	    typedef mpl::integral_c_tag tag;
	    static const int value = A::value + B::value;
	    static const rysq::type max = (A::value >= B::value) ? A_ : B_;
	    static const rysq::type min = (A::value < B::value) ? A_ : B_;
	};

	template<rysq::type A_, rysq::type B_, rysq::type C_, rysq::type D_>
	class braket {
	public:
	    typedef shell<A_> A;
	    typedef shell<B_> B;
	    typedef shell<C_> C;
	    typedef shell<D_> D;

	    typedef state<A_,B_> bra;
	    typedef state<C_,D_> ket;


	    static const int L = bra::L + ket::L;
	    static const int nci = bra::nc0;
	    static const int ncj = bra::nc1;
	    static const int nck = ket::nc0;
	    static const int ncl = ket::nc1;
	    static const int nc = nci*ncj*nck*ncl;
	    static const int size = bra::size*ket::size;

	    typedef mpl::integral_c_tag tag;
	    static const int value = bra::value + ket::value;
	    typedef typename mpl::max<bra,ket>::type max;
	    typedef typename mpl::min<bra,ket>::type min;

	    // 	    typedef util::byteset::type value_type;
	    // 	    static const value_type value = util::byteset::set<A_,B_,C_,D_>::value;
 	};

	template<class A_, class B_, bool apply = true>
	struct transpose {
	    typedef typename mpl::if_<mpl::bool_<apply>, B_, A_>::type A;
	    typedef typename mpl::if_<mpl::bool_<apply>, A_, B_>::type B;
	    static const bool value = apply;
	};

	template<rysq::type A, rysq::type B, bool apply>
	struct transpose <shell<A>, shell<B>, apply> {
	    typedef typename
	    mpl::if_<mpl::bool_<apply>, meta::state<B,A>, meta::state<A,B> >::type type;
	    static const bool value = apply;
	};

	template<rysq::type A_, rysq::type B_, rysq::type C_, rysq::type D_>
	struct sorted {

	    typedef state<A_,B_> bra_;
	    typedef state<C_,D_> ket_;

	    typedef shell<bra_::max> A;
	    typedef shell<bra_::min> B;
	    typedef shell<ket_::max> C;
	    typedef shell<ket_::min> D;

	    typedef meta::braket<A::type, B::type, C::type, D::type> ABCD;
	    typedef meta::braket<C::type, D::type, A::type, B::type> CDAB;

	    static const bool t0_ = shell<A_>::value < shell<B_>::value;
	    static const bool t1_ = shell<C_>::value < shell<D_>::value;
	    typedef mpl::bool_<((A::value == C::value) && (B::value < D::value)) ||
				(A::value < C::value)> t_;
				
	    typedef typename mpl::if_<t_, CDAB, ABCD>::type braket;

	    static const int value = rysq::mpl::transpose<t0_, t1_, t_::value>::value;

	};

	template<class bra, class ket>
	struct sort;

	template<type A, type B, type C, type D>//bra>//, class ket>
	struct sort<state<A,B>, state<C,D> > {
	    typedef typename sorted<A,B,C,D>::braket braket;
	};

    }
}

#endif /* RYSQ_META_HPP */

