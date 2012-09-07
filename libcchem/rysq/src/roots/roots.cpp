/**
   @file
   @brief Rys quadrature roots and weights implementation
*/

#include <cstdlib>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>

#include "roots/roots.hpp"
#include "opq.h"
#include <stdlib.h>
#include <assert.h>

static bool initialized = false;

void rysq::roots_initialize() {
    rysq::detail::stieltjes_init();
#define ROOTS(z, N, text) rysq::asymptotic::initialize<N,double>();
    BOOST_PP_REPEAT_FROM_TO(1, 14, ROOTS, nil);
#undef ROOTS
    initialized = true;
}

void rysq::roots_finalize() {
    rysq::detail::stieltjes_finalize();
    initialized = false;
}

double *rysq::detail::STIELTJES_X[8]; 

double *rysq::detail::STIELTJES_W[8]; 

/**
   @brief Initializes auxiliary Stieltjes quadratures
*/
void rysq::detail::stieltjes_init() {
    double a[55], b[55];

    for(int i = 0; i < 8; ++i) {
	int N = STIELTJES_N[i]; 

	STIELTJES_X[i] = (double*)malloc(N*sizeof(double));
	STIELTJES_W[i] = (double*)malloc(N*sizeof(double));

	a[0] = 0.5;
	b[0] = 1.0; 

	for(int j = 1; j < N; ++j) {
	    a[j] = 0.5;
	    b[j] = 0.25/(4.0 - (1.0/(j*j)));
	}

	opq::coefficients(N, a, b, STIELTJES_X[i], STIELTJES_W[i]); 
    }
}

/**
   @brief Finalizes auliliary Stieltjes quadratures
*/
void rysq::detail::stieltjes_finalize() {
    
    for(int i = 0; i < 8; ++i) {

	if(STIELTJES_X[i]) {
	    free(STIELTJES_X[i]);
	    STIELTJES_X[i] = NULL;
	}

	if(STIELTJES_W[i]) {
	    free(STIELTJES_W[i]);
	    STIELTJES_W[i] = NULL;
	}

    }
}

void rysq_roots_initialize() {
    rysq::roots_initialize();
}


#ifdef INTEGER64
void rysq_roots(int64_t *N, double *X, double *t, double *W)
#else
void rysq_roots(int32_t *N, double *X, double *t, double *W)
#endif
{    
#define ROOTS(z, N_, text) if (*N == N_) {	\
	double t_[N_], W_[N_+(N_==0)];		\
	rysq::roots<N_>(*X, t_, W_);		\
	for (int i = 0; i < N_; ++i) {		\
	    t[i] = t_[i];			\
	    W[i] = W_[i];			\
	}					\
	if (N_ == 0) *W = W_[0];		\
    }						\
    else

    BOOST_PP_REPEAT_FROM_TO(0, 14, ROOTS, nil) {
	fprintf(stderr, "rysq::roots: too many roots %i\n", int(*N));
    }

#undef ROOTS

}

#pragma weak rysq_roots_initialize_ = rysq_roots_initialize
#pragma weak rysq_roots_ = rysq_roots

