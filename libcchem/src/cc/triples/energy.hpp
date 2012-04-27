#ifndef CC_TRIPLES_ENERGY_HPP
#define CC_TRIPLES_ENERGY_HPP

#include "cc/cc.hpp"

#include <utility>
#include <boost/typeof/typeof.hpp>
#include <boost/numeric/ublas/adaptor.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/multi_array/multi_array_ref.hpp>
#include <boost/array.hpp>

#include <boost/fusion/include/pair.hpp>
#include <boost/fusion/include/map.hpp>
#include <boost/fusion/include/intrinsic.hpp>
#include <boost/fusion/include/for_each.hpp>

#include "foreach.hpp"
#include "array/hdf5.hpp"
#include "array/permute.hpp"

#include "blas.hpp"

namespace cc {
namespace triples {
namespace detail {

    struct Energy {    	

	template<int I>
	struct Index {
	    Index(int i) : value_(i) {}
	    operator int&() { return value_; }
	    operator const int&() const { return value_; }
	private: int value_;
	};

	struct Range {
	    const int start, stop, size;
	    Range(int start, int stop)
		: start(start), stop(stop), size(stop - start) {}
	};

	template<int N>
	struct Tile {
	    typedef double A1[N];
	    typedef double A2[N][N];
	    typedef double A3[N][N][N];
	    A1 ai, bj, ck;
	    A2 abij, acik, bcjk;
	    double cba[N][N][N];
	    double cab[N][N][N];
	    double bca[N][N][N];
	    double bac[N][N][N];
	    double abc[N][N][N];
	    double acb[N][N][N];

	    template<class A>
	    struct Tie {
		Tie(A a, A b, A c) : _0(a), _1(b), _2(c) {}
		template<int I>
		A get() const { return get(int_<I>()); }
	    private:
		template<int I> struct int_ {};
		A get(int_<0>) const { return _0; }
		A get(int_<1>) const { return _1; }
		A get(int_<2>) const { return _2; }
		A _0, _1, _2;
	    };

	    template<class T1, class T2, class T3>
	    Tile(const Range &ra, const Range &rb, const Range &rc,
		 int i, int j, int k,
		 const T1 &t1, const T2 &t2, const T3 &t3) {
		load(ra, i, t1, ai);
		load(rb, j, t1, bj);
		load(rc, k, t1, ck);

		load(ra, rb, i, j, t2.ij, abij);
		load(ra, rc, i, k, t2.ik, acik);
		load(rb, ra, j, k, t2.jk, bcjk);

		Tie<const Range&> r(ra, rb, rc);
		load<0,1,2>(r, t3, abc);
		load<0,2,1>(r, t3, acb);
		load<1,0,2>(r, t3, bac);
		load<1,2,0>(r, t3, bca);
		load<2,0,1>(r, t3, cab);
		load<2,1,0>(r, t3, cba);
	    }
		
	    template<class T>
	    static void load(const Range& r, int i, const T &t, A1 &A) {
		for (int a = 0; a < r.size; ++a) {
		    A[a] = t(a,i);
		}
	    }

	    template<class T>
	    static void load(const Range& ra, const Range& rb,
			     int i, int j, const T &t, A2 &A) {
		for (int b = 0; b < rb.size; ++b) {
		    for (int a = 0; a < ra.size; ++a) {
			A[b][a] = t(a,b);
		    }
		}
	    }

	    template<int I, int J, int K, class T3>
	    static void load(Tie<const Range&> r, const T3 &t3, A3 &A) {
		const Range& ra = r.template get<I>();
		const Range& rb = r.template get<J>();
		const Range& rc = r.template get<K>();
		const int I_ = (I == 0) ? 0 : ((J == 0) ? 1 : 2);
		const int J_ = (I == 1) ? 0 : ((J == 1) ? 1 : 2);
		const int K_ = (I == 2) ? 0 : ((J == 2) ? 1 : 2);
		for (int c = 0; c < rc.size; ++c) {
		    for (int b = 0; b < rb.size; ++b) {
			BOOST_AUTO(t3_a,
				   t3[c+rc.start][b+rb.start].begin()+ra.start);
			for (int a = 0; a < ra.size; ++a) {
			    Tie<int&> i(a,b,c);
			    A[i.template get<K_>()]
				[i.template get<J_>()]
				[i.template get<I_>()] = *t3_a++;
			    
			}
		    }
		}
	    }

	};

	template<int N, class T3>
	static
	Correction evaluate(const Matrix &t1,
			    Permutations<const Matrix&,3> t2,
			    Permutations<const Matrix&,6> vvoo,
			    const T3 &t3,
			    const Vector &eh, const Vector &ep,
			    Matrix &u1,
			    Index<'i'> i, Index<'j'> j, Index<'k'> k,
			    Range ra, Range rb, Range rc) {

	    const double symmetry = ((i == j || j == k) ? 0.5 : 1);
	    double dijk = eh(i) + eh(j) + eh(k);

#define t1(a,i) t. a ## i [a]
#define t2(a,b,i,j) (t. a ## b ## i ## j [b][a])
#define t3(a,b,c) t. a ## b ## c [c_][b_][a_]
#define V(a,b,i,j) (vvoo. i ## j (a + r ## a.start, b + r ## b.start))

	    Tile<N> t(ra, rb, rc, i, j, k, t1, t2, t3);

	    Correction C;
	    for (int c = 0; c < rc.size; ++c) {
		for (int b = 0; b < rb.size; ++b) {
		    double dbc = ep(b+rb.start) + ep(c+rc.start);
		    double t1_jb = t1(b,j);
		    double t1_kc = t1(c,k);
		    double t1_jbkc = t1_jb*t1_kc;

		    double Vij = V(b,c,i,j);
		    double Vik = V(b,c,i,k);
		    double Vjk = V(b,c,j,k);
		    double Vji = V(b,c,j,i);
		    double Vki = V(b,c,k,i);
		    double Vkj = V(b,c,k,j);

		    for (int a = 0; a < ra.size; ++a) {
			int &a_ = a;
			int &b_ = b;
			int &c_ = c;

			double s = symmetry*(!(a+ra.start == b+rb.start &&
					       b+rb.start == c+rc.start));

			double d = s/(dijk - (ep(a+ra.start) + dbc));
			double abc1 = t1(a,i)*t1_jbkc;
			double abc2 = (t1(a,i)*t2(b,c,j,k) +
				       t1_jb*t2(a,c,i,k) +
				       t1_kc*t2(a,b,i,j));
			double abc3 = t3(a,b,c);
			double f = (8*abc3 -
				    4*(t3(a,c,b)+t3(c,b,a)+t3(b,a,c)) +
				    2*(t3(b,c,a)+t3(c,a,b)))*d;
			C.ots += f*abc1;
			C.otd += f*abc2;
			C.etd += f*abc3;

			double t3_[] =  {
			    2*(t3(a,b,c) - t3(b,a,c)) - t3(a,c,b) + t3(b,c,a),
			    2*(t3(a,c,b) - t3(b,c,a)) - t3(a,b,c) + t3(b,a,c),
			    2*(t3(b,a,c) - t3(a,b,c)) - t3(c,a,b) + t3(c,b,a),
			    2*(t3(b,c,a) - t3(a,c,b)) - t3(c,b,a) + t3(c,a,b),
			    2*(t3(c,a,b) - t3(c,b,a)) - t3(b,a,c) + t3(a,b,c),
			    2*(t3(c,b,a) - t3(c,a,b)) - t3(b,c,a) + t3(a,c,b)
			};

			// ET[S]
			u1(a,i) += d*(t3_[0]*Vjk + t3_[1]*Vkj);
			u1(a,j) += d*(t3_[2]*Vik + t3_[4]*Vki);
			u1(a,k) += d*(t3_[3]*Vij + t3_[5]*Vji);
		    }
		}
	    }

#undef t1
#undef t2
#undef t3
#undef V
	    return C;
	}
    };

} // namespace detail
} // namespace triples
} // namespace cc

#endif // CC_TRIPLES_ENERGY_HPP
