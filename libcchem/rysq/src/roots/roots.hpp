#ifndef RYSQ_ROOTS_HPP
#define RYSQ_ROOTS_HPP


#include <boost/config.hpp>

#ifdef __CUDACC__
#define constant__ __device__ const
#else
#define constant__ static const
#endif

#ifndef BOOST_GPU_ENABLED
#define BOOST_GPU_ENABLED
#endif

#ifndef __CUDACC__
#include "roots/asymptotic.hpp"
#include "roots/opq.h"
#endif

#include <boost/utility/enable_if.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/at.hpp>

#include <stdio.h>
#include <math.h>
#include <stdint.h>

namespace rysq {
namespace detail {

    struct serial_tag {};

    template<class T, class R = void>
    struct disable_if_serial {
	typedef R type;
    };

    template<class R>
    struct disable_if_serial<serial_tag, R> {}; 

#include "roots/evaluate.hpp"

    template<size_t N, typename T>
    struct roots;

#define TYPE__ float
#include "roots/roots_type.hpp"
#undef TYPE__

#define TYPE__ double
#include "roots/roots_type.hpp"
#undef TYPE__

#ifndef __CUDACC__

    typedef long double long_double;
#define TYPE__ long_double
#include "roots/roots_type.hpp"
#undef TYPE__

#endif

}
}

namespace rysq {

    void roots_initialize();
    void roots_finalize();

    template<size_t N, typename T>
    BOOST_GPU_ENABLED
    typename boost::enable_if_c<(N == 0)>::type
    roots(T X, T (&t2)[0], T (&W)[1]) {
	return detail::roots<0,T>::evaluate(X, t2, W);
    }

    template<size_t N, typename T>
    BOOST_GPU_ENABLED
    void roots(T X, T (&t2)[N], T (&W)[N]) {
	return detail::roots<N,T>::evaluate(X, t2, W);
    }

    template<size_t N, typename T>
    BOOST_GPU_ENABLED
    void roots(T X, T *t2, T *W, const unsigned short thread) {
	detail::roots<N,T>::evaluate(X, t2, W, thread);
    }

}

#ifndef __CUDACC__

namespace rysq {
namespace detail {

    /**
       @brief Auxiliary quadrature size corresponding to a number of roots
    */
    template<size_t n>
    struct STIELTJES {
    	typedef typename boost::mpl::vector_c
    	<int,0, 1, 2, 2, 3, 4, 4, 4, 5, 6, 6, 7, 7>::type index_vector;
    	static const size_t index = boost::mpl::at_c<index_vector,n-1>::type::value;
    	static const size_t N = 20 + index*5;
    };

    static const int STIELTJES_N[] = { 20, 25, 30, 35, 40, 45, 50, 55 }; 

    /**
       @brief Auxiliary quadrature index corresponding to a number of roots
    */
    static const int STIELTJES_INDEX[] = { 0, 1, 2, 2, 3, 4, 4, 4, 5, 6, 6, 7, 7 };

    /**
       @brief Auxiliary quadrature roots
    */ 
    extern double *STIELTJES_X[8]; 

    /**
       @brief Auxiliary quadrature weights
    */
    extern double *STIELTJES_W[8]; 

    void stieltjes_init();

    void stieltjes_finalize();

    /**
       @brief Evaluates Rys quadrature roots and weights using Stieltjes method
       @param n Number of roots
       @param X X value
       @param[out] x Roots
       @param[out] w Weights 
    */   
    template<int n, typename T>
    int stieltjes(T X, T *x, T *w) {
    
	// int is = STIELTJES_INDEX[n-1];
	// int N = STIELTJES_N[is];
	static const int N = STIELTJES<n>::N;
	static const int is = STIELTJES<n>::index;

	T rg[N], wg[N];
	T a[n], b[n];
	
	for(int i = 0; i < N; ++i) {
	    T S = STIELTJES_X[is][i];
	    rg[i] = S*S;
	    wg[i] = STIELTJES_W[is][i]*exp(-X*S*S);
	}

	int status = 0;
	status = opq::stieltjes<n,N>(rg, wg, a, b);
	if (status != 0) {
	    fprintf(stderr, "rysq::roots:  opq::stieltjes returned %i\n", status);
	}

	status = opq::coefficients(n, a, b, x, w);
	if (status != 0) {
	    fprintf(stderr, "rysq::roots:  opq::coefficients returned %i\n", status);
	}
		     
	//    for(int i = 0; i < n; ++i) {
	//	x[i] = x[i]/(1.0 - x[i]);
	// }

	return status;
 
    }

    template<size_t N, typename T>
    struct roots {
	BOOST_STATIC_ASSERT((N > 0));
	static void evaluate(T X, T *t2, T *W) {
	    const double asymp[] = 
		{ 29, 37, 43, 49, 55, 60, 65, 71, 76, 81,
		  86, 91, 96, 100, 100, 100, 100, 100, 100, 100 };
	    if (X > asymp[(N ? N-1 : 0)])
		asymptotic::roots<N>(X, t2, W);
	    else stieltjes<N>(X, t2, W);
	}
    };


} // namespace roots
} // namespace rysq

#endif // __CUDACC__

#ifndef __CUDACC__
#undef constant__
#endif

extern "C" {

    void rysq_roots_initialize();
#ifdef INTEGER64
    void rysq_roots(int64_t *N, double *X, double *t, double *W);
#else
    void rysq_roots(int32_t *N, double *X, double *t, double *W);
#endif

    void rysq_roots_initialize_();
#ifdef INTEGER64
    void rysq_roots_(int64_t *N, double *X, double *t, double *W);
#else
    void rysq_roots_(int32_t *N, double *X, double *t, double *W);
#endif

}

#endif //  RYSQ_ROOTS_HPP
