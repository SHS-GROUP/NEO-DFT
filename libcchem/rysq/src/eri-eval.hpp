#ifndef _RYSQ_ERI_EVAL_HPP_
#define _RYSQ_ERI_EVAL_HPP_

#include "rysq/core.hpp"
#include "kernels/quadrature.h"
#include "kernels/functor.hpp"
#include "vector.hpp"
#include "roots/rysq_roots.h"
#include "transpose.hpp"

namespace rysq {
    namespace eri {

	template <class braket>
	static void eval(const Quartet<Shell> &quartet,
			 const Vector<3> &ri, const Vector<3> &rj,
			 const Vector<3> &rk, const Vector<3> &rl,
			 double *I, double scale, double cutoff);

    }
}



template<class braket>
void rysq::eri::eval(const Quartet<Shell> &quartet,
		     const Vector<3> &ri, const Vector<3> &rj,
		     const Vector<3> &rk, const Vector<3> &rl,
		     double *I, double scale, double cutoff) {

    typedef kernels::braket<braket::A::type, braket::B::type, 
	braket::C::type, braket::D::type> kernel;

    const Shell &a = quartet[0];
    const Shell &b = quartet[1];
    const Shell &c = quartet[2];
    const Shell &d = quartet[3];

    static const int L = braket::L;
    static const int N = (L > 0) ? L/2 + 1 : 0;

    Vector<3> rij = ri - rj;
    Vector<3> rkl = rk - rl;

    double rij2 = rij.inner();
    double rkl2 = rkl.inner();

    double eij[a.K*b.K] __attribute__ ((aligned(16)));
    for(int Kj = 0, Kij = 0; Kj < b.K; ++Kj) {
	for(int Ki = 0; Ki < a.K; ++Ki, ++Kij) {
	    double ai = a(Ki);
	    double aj = b(Kj);
	    eij[Kij] = exp(-ai*aj*rij2/(ai + aj));
	}
    }

    double Q[braket::size]  __attribute__ ((aligned(16))) = { 0.0 };

    for (int Kl = 0; Kl < d.K; ++Kl) {
	for (int Kk = 0; Kk < c.K; ++Kk) {
	    double ak = c(Kk);
	    double al = d(Kl);
	    double B = ak + al;

	    Vector<3> rB = Vector<3>::center(ak, rk, al, rl);
	    Vector<3> rBk = rB - rk;

	    double ekl = exp(-ak*al*rkl2/B);

	    static const int nckl = braket::nck*braket::ncl;
	    double Ckl[nckl] __attribute__ ((aligned(16)));

	    for(int l = 0; l < braket::ncl; ++l) {
		for(int k = 0; k < braket::nck; ++k) {
		    double q = c(Kk,k)*d(Kl,l);
		    Ckl[k+l*braket::nck] = q;
		}
	    }
		    
	    for(int Kj = 0; Kj < b.K; ++Kj) {
		for(int Ki = 0; Ki < a.K; ++Ki) {
		    double ai = a(Ki);
		    double aj = b(Kj); 
		    double A = ai + aj;
		    double e = eij[Ki+Kj*a.K]*ekl;
		    double eAB = e/(A*B*sqrt(A+B));

		    double C[braket::nc] __attribute__ ((aligned(16)));
		    double Cmax = 0.0;

		    for(int kl = 0, ijkl = 0; kl < nckl; ++kl) {
			double q = eAB*Ckl[kl];
			for(int j = 0; j < braket::ncj; ++j) {
			    for(int i = 0; i < braket::nci; ++i, ++ijkl) {
				double p = q*a(Ki,i)*b(Kj,j);
				C[ijkl] = p;
				Cmax = std::max(fabs(p), Cmax);
			    }
			}
		    }

		    if (Cmax < cutoff) continue;

		    Vector<3> rA = Vector<3>::center(ai, ri, aj, rj);
		    Vector<3> rAi = rA - ri;
		    Vector<3> rAB = rA - rB;

		    double rho = (A*B)/(A + B);
		    double X =  rho*rAB.inner();
		    Vector<N> t2;
		    Vector<N+(L == 0)> W;
		    rysq::roots<N>(X, t2, W);

		    t2 /= (A + B);
		    kernel::eval(A, B, rAi, rBk, rAB, rij, rkl, t2, W, C, Q);

		}
	    }
	}
    }
    kernel::reorder(Q);

    if (quartet.transpose_mask) {
	Transpose::apply<braket::A::size, braket::B::size,
	    braket::C::size, braket::D::size>(Q, quartet.transpose_mask);
    }
    util::scale(braket::size, I, scale, Q);

}


#endif /* _RYSQ_ERI_EVAL_HPP_ */
