#include "opq.h"

#include <assert.h>
#include <stdio.h>
#include <boost/math/constants/constants.hpp>

namespace rysq {
namespace asymptotic {

    template<size_t N, typename T>
    struct asymptotic_ {
	static void initialize() {
	    static const size_t n = 2*N;
	    T beta[n], alpha[n] = { 0 };
	    beta[0] = boost::math::constants::root_pi<T>();
	    for (size_t i = 1; i < n; ++i) {
		beta[i] = T(i)/2;
	    }
	    T r[n], w[n];
	    int status = opq::coefficients(n, alpha, beta, r, w);
	    if (status != 0)  {
		fprintf(stderr, "%s: opq::coefficients returned %i\n",
			__FILE__, status);
		return;
	    }
	    // CALL RYSGW_(N,ALPHA,BETA,EPS,RTS,W_TS,IERR,W_RK) for
	    for (size_t i = 0; i < N; ++i) {
		size_t j = i + N;
		R_[i] = r[j]*r[j];
		W_[i] = w[j];
	    }
	    initialized_ = true;
	}
	static void roots(T X, T *R, T *W) {
	    assert(initialized_);
	    T r = 1/X;
	    T w = 1/sqrt(X);
	    for (size_t i = 0; i < N; ++i) {
		R[i] = r*R_[i];
		W[i] = w*W_[i];
	    }
	}
    private:
	static T R_[N], W_[N];
	static bool initialized_;
    };
    template<size_t N, typename T> T asymptotic_<N,T>::R_[N];
    template<size_t N, typename T> T asymptotic_<N,T>::W_[N];
    template<size_t N, typename T> bool asymptotic_<N,T>::initialized_;

    template<size_t N, typename T>
    void roots(T X, T *R, T *W) {
	asymptotic_<N,T>::roots(X, R, W);
    }

    template<size_t N, typename T>
    void initialize() {
	asymptotic_<N,T>::initialize();
    }

} // namespace asymptotic
} // namespace rysq

