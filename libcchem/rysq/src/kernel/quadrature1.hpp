#ifndef RYSQ_KERNEL_QUADRATURE1_HPP
#define RYSQ_KERNEL_QUADRATURE1_HPP

#include "kernel/quadrature.hpp"
#include "kernel/quadrature1-impl.hpp"

#include <math.h>
#include "rysq-core.hpp"
#include "vector.hpp"

#include "roots/roots.hpp"
#include "kernel/transfer.hpp"
#include "kernel/recurrence.hpp"


namespace rysq {
namespace kernel {
namespace quadrature {


template<class bra, size_t N, size_t NT,
	 class Transform, typename T, typename U>
static void apply(rysq::type c, rysq::type d,
		  const Primitives<T,U> &primitives,
		  Transform &transform, double scale);

template<class bra, size_t N,template<typename, size_t> class align,
	 class Transform, typename T, typename U>
static void apply(const Quartet<Shell> &quartet,
		  const vector<3> &ri, const vector<3> &rj,
		  const vector<3> &rk, const vector<3> &rl,
		  double scale, double cutoff,
		  Primitives<T,U> &primitives,
		  Transform &transform) {

    static const size_t NT = N + align<T,N>::value;
    //static const size_t NU = N + align<U,N>::value;

    const Shell &a = quartet[0];
    const Shell &b = quartet[1];
    const Shell &c = quartet[2];
    const Shell &d = quartet[3];

    vector<3> rij = ri - rj;
    vector<3> rkl = rk - rl;

    double rkl2 = rkl.dot();
    double rij2 = rij.dot();

    double eij[a.K*b.K];
    for (int Kj = 0, Kij = 0; Kj < b.K; ++Kj) {
	for (int Ki = 0; Ki < a.K; ++Ki, ++Kij) {
	    double A = a(Ki) + b(Kj);
	    double A1 = 1.0/A;
	    eij[Kij] = exp(-a(Ki)*b(Kj)*A1*rij2);
	}
    }

    int K = 0;

    for (int Kl = 0; Kl < d.K; ++Kl) {
	for (int Kk = 0; Kk < c.K; ++Kk) {
	    double ak = c(Kk);
	    double al = d(Kl);
	    double B = ak + al;

	    vector<3> rB = vector<3>::center(ak, rk, al, rl);
	    vector<3> rBk = rB - rk;

	    double ekl = exp(-ak*al*rkl2/B);

	    double Ckl[4] __attribute__ ((aligned(16)));
	    double Ckl_max = 0.0;
	    for(int l = 0, kl = 0; l < d.nc; ++l) {
		for(int k = 0; k < c.nc; ++k, ++kl) {
		    Ckl[kl] = c(Kk,k)*d(Kl,l);
		    Ckl_max = std::max(Ckl_max, fabs(Ckl[kl]));
		}
	    }
		    
	    for(int Kj = 0, Kij = 0; Kj < b.K; ++Kj) {
		for(int Ki = 0; Ki < a.K; ++Ki, ++Kij) {
		    double ai = a(Ki);
		    double aj = b(Kj); 
		    double A = ai + aj;
		    double e = eij[Ki+Kj*a.K]*ekl;
		    double eAB = e/(A*B*sqrt(A+B));

		    double Cij[bra::nc] __attribute__ ((aligned(16)));
		    double Cij_max = 0.0;

		    for(int j = 0, ij = 0; j < bra::B::nc; ++j) {
			for(int i = 0; i < bra::A::nc; ++i, ++ij) {
			    Cij[ij] = eAB*a(Ki,i)*b(Kj,j);
			    Cij_max = std::max(fabs(Cij[ij]), Cij_max);
			}
		    }
		    if (Cij_max*Ckl_max < cutoff) continue;

		    for(int kl = 0; kl < (c.nc*d.nc); ++kl) {
			for(int ij = 0; ij < bra::nc; ++ij) {
			    primitives.C[kl][ij + K*bra::nc] = Cij[ij]*Ckl[kl];
			}
		    }

		    vector<3> rA = vector<3>::center(ai, ri, aj, rj);
		    vector<3> rAi = rA - ri;
		    vector<3> rAB = rA - rB;

		    double rho = A*B/(A + B);
		    U X =  rho*rAB.dot();
		    vector<N,U> W, t2;
		    rysq::roots<N>(X, t2.elems, W.elems);
		    t2 /= (A + B);

		    static const int LA = bra::A::L;
		    static const int LB = bra::B::L;

		    T *Ix = primitives.Ix(K);
		    T *Iy = primitives.Iy(K);
		    T *Iz = primitives.Iz(K);

		    typedef kernel::recurrence<bra::L,N,align> recurrence;
		    typedef kernel::transfer<LA,LB,N,align> transfer;

		    if (LB + d.L == 0) {
		    	recurrence::apply(c.L + d.L, A, B, rAB, rAi, rBk, 
					  t2.data(), W.data(), Ix, Iy, Iz);
		    }
		    else {
		    	U *Gx = primitives.Gx(K);
		    	U *Gy = primitives.Gy(K);
		    	U *Gz = primitives.Gz(K);
		    	U *tmp = primitives.transfer();

			recurrence::apply(c.L + d.L, A, B, rAB, rAi, rBk, 
		    	 		  t2.data(), W.data(), Gx, Gy, Gz);

		    	transfer::apply(c.L, d.L, U(rij[0]), U(rkl[0]), Gx, Ix, tmp);
		     	transfer::apply(c.L, d.L, U(rij[1]), U(rkl[1]), Gy, Iy, tmp);
		     	transfer::apply(c.L, d.L, U(rij[2]), U(rkl[2]), Gz, Iz, tmp);
		    }

		    ++K;

		    // contract primitives
		    if (K == primitives.Kmax) {
			primitives.K = K;
			apply<bra,N,NT>(c, d, primitives, transform, scale);
			K = 0;
		    }

		}
	    }
	}
    }

    // Final contraction
    primitives.K = K;
    if (primitives.K) apply<bra,N,NT>(c, d, primitives, transform, scale);

}

template<class bra, size_t N, size_t NT,
	 class Transform, typename T, typename U>
static void apply(rysq::type c, rysq::type d,
		  const Primitives<T,U> &primitives,
		  Transform &transform, double scale) {

    static const int Nij = NT*(bra::A::L+1)*(bra::B::L+1);

    const int Lc = abs(c);
    const int Ld = abs(d);

    const int dim2d = Nij*(Lc+1)*(Ld+1);

    const int K = primitives.K;
    const double* const *C = primitives.C;
    const T *Ix = primitives.Ix(0);
    const T *Iy = primitives.Iy(0);
    const T *Iz = primitives.Iz(0);

    int spk = (c < 0);
    int spl = (d < 0);

    const int c_first = shell(c).begin();
    const int d_first = shell(d).begin();
    const int c_last = shell(c).end() - 1;
    const int d_last = shell(d).end() - 1;

    for(int l = d_first, kl = 0; l <= d_last; ++l) {
	const int lx = (Lc+1)*LX[l];
	const int ly = (Lc+1)*LY[l];
	const int lz = (Lc+1)*LZ[l];

	const int lsp = (spl && l) << spk;

	for(int k = c_first; k <= c_last; ++k, ++kl) {
	    const double *Ckl = C[(spk && k) + lsp];

	    const int klx = Nij*(lx + LX[k]);
	    const int kly = Nij*(ly + LY[k]);
	    const int klz = Nij*(lz + LZ[k]);

	    static const int ni = bra::A::size;
	    static const int nj = bra::B::size;
	    double I[ni*nj] __attribute__((aligned(16))) = { 0.0 };
	    
	    int flags = 0;
	    double screen = 0.0;

	     // std::cout << size_t(&Ix[klx])%16
	     // 	      << size_t(&Iy[kly])%16
	     // 	       << size_t(&Iz[klz])%16 << NT << std::endl;

	    typedef quadrature::impl<bra> impl;
	    impl::template apply<N,T,NT>(flags, screen, K, Ckl, dim2d,
	     				 &Ix[klx], &Iy[kly], &Iz[klz], 1.0, I);

	    transform(k - c_first, l - d_first, kl, I, scale);
	}
    }
}

} // namespace rysq
} // namespace kernel
} // namespace quadrature

#endif // RYSQ_KERNEL_QUADRATURE1_HPP

