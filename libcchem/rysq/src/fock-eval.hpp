#ifndef RYSQ_FOCK_EVAL_HPP
#define RYSQ_FOCK_EVAL_HPP

#include "rysq-fock.hpp"
#include "fock-transform.hpp"
#include "meta.hpp"

namespace rysq {
namespace fock {

typedef Fock::density_set density_set;
typedef Fock::fock_set fock_set;

template<class braket>
inline void eval(const density_set &D, const fock_set &F,
		 const double *__restrict Eri,
		 const Fock::Parameters &parameters) {

    // using util::add;
    // using util::copy;

    double xscale = parameters.xscale;
    double cscale = parameters.cscale;

    static const int ni = braket::A::size;
    static const int nj = braket::B::size;
    static const int nk = braket::C::size;
    static const int nl = braket::D::size;


    double Dij[ni*nj] __attribute__ ((aligned(16)));
    double Dkl[nk*nl] __attribute__ ((aligned(16)));
    double Dik[ni*nk] __attribute__ ((aligned(16)));
    double Dil[ni*nl] __attribute__ ((aligned(16)));
    double Djk[nj*nk] __attribute__ ((aligned(16)));
    double Djl[nj*nl] __attribute__ ((aligned(16)));

    copy(ni*nj, Dij, D[0]);
    copy(nk*nl, Dkl, D[1]);
    copy(ni*nk, Dik, D[2]);
    copy(ni*nl, Dil, D[3]);
    copy(nj*nk, Djk, D[4]);
    copy(nj*nl, Djl, D[5]);

    double Fij[ni*nj] __attribute__ ((aligned(16))) = { 0.0 };
    double Fkl[nk*nl] __attribute__ ((aligned(16))) = { 0.0 };
    double Fik[ni*nk] __attribute__ ((aligned(16))) = { 0.0 };
    double Fil[ni*nl] __attribute__ ((aligned(16))) = { 0.0 };
    double Fjk[nj*nk] __attribute__ ((aligned(16))) = { 0.0 };
    double Fjl[nj*nl] __attribute__ ((aligned(16))) = { 0.0 };

    for (int l = 0, kl = 0, ijkl = 0; l < nl; ++l) {

	for (int k = 0; k < nk; ++k, ++kl) {
	    int jk = k*nj;
	    int jl = l*nj;

	    for (int j = 0, ij = 0; j < nj; ++j, ++jk, ++jl) {
		int ik = k*ni;
		int il = l*ni;

		for (int i = 0; i < ni; ++i, ++ij, ++ik, ++il, ++ijkl) {
		    //double q = scale*Eri[ijkl];
		    double q = Eri[ijkl];
		    //if(fabs(q) < 1e-11) continue;

		    Fij[ij] += Dkl[kl]*q;
		    Fkl[kl] += Dij[ij]*q;
		    Fik[ik] += Djl[jl]*q;
		    Fil[il] += Djk[jk]*q;
		    Fjk[jk] += Dil[il]*q;
		    Fjl[jl] += Dik[ik]*q;
		}
	    }

	}
    }
    add(ni*nj, F[0], cscale, Fij);
    add(nk*nl, F[1], cscale, Fkl);
    add(ni*nk, F[2], xscale, Fik);
    add(ni*nl, F[3], xscale, Fil);
    add(nj*nk, F[4], xscale, Fjk);
    add(nj*nl, F[5], xscale, Fjl);
}

template<int ni, int nj, int nk, int nl>
inline void eval(const Data::density_type &D, const Data::fock_type &F,
		 const double *__restrict Q, double scale) {

    const double * __restrict Dij = D[0];
    const double * __restrict Dkl = D[1];
    const double * __restrict Dik = D[2];
    const double * __restrict Dil = D[3];
    const double * __restrict Djk = D[4];
    const double * __restrict Djl = D[5];

    double * __restrict Fij = F[0];
    double * __restrict Fkl = F[1];
    double * __restrict Fik = F[2];
    double * __restrict Fil = F[3];
    double * __restrict Fjk = F[4];
    double * __restrict Fjl = F[5];

    for (int l = 0, kl = 0, ijkl = 0; l < nl; ++l) {

	for (int k = 0; k < nk; ++k, ++kl) {
	    int jk = k*nj;
	    int jl = l*nj;

	    for (int j = 0, ij = 0; j < nj; ++j, ++jk, ++jl) {
		int ik = k*ni;
		int il = l*ni;

		for (int i = 0; i < ni; ++i, ++ij, ++ik, ++il, ++ijkl) {
		    double q = scale*Q[ijkl];
		    // std::cout << Dkl[kl] << std::endl;
		    Fij[ij] += Dkl[kl]*q;
		    Fkl[kl] += Dij[ij]*q;
		    Fik[ik] += Djl[jl]*q;
		    Fil[il] += Djk[jk]*q;
		    Fjk[jk] += Dil[il]*q;
		    Fjl[jl] += Dik[ik]*q;
		}
	    }

	}
    }
}

template<int M, int N>
void eval(int k, int l, int kl,
	  const Data::density_type &D, const Data::fock_type &F,
	  const double * __restrict I, double scale) {

    const double * __restrict Dij = D[0];
    const double * __restrict Dkl = D[1];
    const double * __restrict Dik = D[2];
    const double * __restrict Dil = D[3];
    const double * __restrict Djk = D[4];
    const double * __restrict Djl = D[5];

    double * __restrict Fij = F[0];
    double * __restrict Fkl = F[1];
    double * __restrict Fik = F[2];
    double * __restrict Fil = F[3];
    double * __restrict Fjk = F[4];
    double * __restrict Fjl = F[5];

    int jk = k*N;
    int jl = l*N;
    for (int j = 0, ij = 0; j < N; ++j, ++jk, ++jl) {

    	int ik = k*M;
    	int il = l*M;
    	for (int i = 0; i < M; ++i, ++ij, ++ik, ++il) {
    	    double q = scale*I[ij];
    	    Fij[ij] += Dkl[kl]*q;
    	    Fkl[kl] += Dij[ij]*q;
    	    Fik[ik] += Djl[jl]*q;
    	    Fil[il] += Djk[jk]*q;
    	    Fjk[jk] += Dil[il]*q;
	    Fjl[jl] += Dik[ik]*q;
    	}
    }

}


namespace {

    template<class P>
    struct Scale;

    template<>
    struct Scale<Fock::Parameters> {
	explicit Scale(const Fock::Parameters &p)
	    : c_(p.cscale), x_(p.xscale) {}
	template<typename T>
	T x(const T &q) const { return x_*q; }
	template<typename T>
	T c(const T &q) const { return c_*q; }
	template<typename T>
	const T& operator()(const T &q) const { return q; }
    private:
	double c_, x_;
    };

    template<rysq::type A, rysq::type B, class S>
    void eval(const rysq::type &c, const rysq::type &d,
	      const density_set &D, const fock_set &F,
	      const double *eri, const S &scale) {

	static const int ni = rysq::meta::shell<A>::size;
	static const int nj = rysq::meta::shell<B>::size;
	const int nk = rysq::shell(c).size;
	const int nl = rysq::shell(d).size;

	for (int l = 0, kl = 0, ijkl = 0; l < nl; ++l) {

	    for (int k = 0; k < nk; ++k, ++kl) {
		int jk = k*nj;
		int jl = l*nj;

		for (int j = 0, ij = 0; j < nj; ++j, ++jk, ++jl) {
		    int ik = k*ni;
		    int il = l*ni;

		    for (int i = 0; i < ni; ++i, ++ij, ++ik, ++il, ++ijkl) {
			double q = scale(eri[ijkl]);
			F[0][ij] += scale.c(D[1][kl]*q);
			F[1][kl] += scale.c(D[0][ij]*q);
			F[2][ik] += scale.x(D[5][jl]*q);
			F[3][il] += scale.x(D[4][jk]*q);
			F[4][jk] += scale.x(D[3][il]*q);
			F[5][jl] += scale.x(D[2][ik]*q);
		    }
		}

	    }
	}
    }

} // namespace

template<rysq::type A, rysq::type B>
void eval(const rysq::type &c, const rysq::type &d,
	  const density_set &D, const fock_set &F,
	  const double *eri, const Fock::Parameters &p) {
    eval<A,B>(c, d, D, F, eri, Scale<Fock::Parameters>(p));
}

} // namespace fock
} // namespace rysq

#endif /* RYSQ_FOCK_EVAL_HPP */
