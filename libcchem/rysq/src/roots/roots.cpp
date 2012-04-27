/**
   @file
   @brief Rys quadrature roots and weights implementation
*/

#include <cstdlib>

#include "roots/roots.hpp"
#include "opq.h"
#include <stdlib.h>
#include <assert.h>
#include <stdexcept>

static bool initialized = false;

void rysq::roots_initialize() {
    rysq::detail::stieltjes_init();
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

