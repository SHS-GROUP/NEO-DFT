#ifndef _UTIL_H
#define _UTIL_H

#include <algorithm>
#include <iostream>
#include <stdint.h>
#include <boost/array.hpp>

#ifdef __CUDACC__
#define _rysq_host_device_ __host__ __device__
#else
#define _rysq_host_device_
#endif

/**
   @file
   @brief Various utility functions
   @internal
*/

/**
   @brief Raises a number to the second power
*/
//#define pow2(x) ((x)*(x))

/**
   @brief Prints error message msg and exits with -1
*/
#define error(msg) {printf("librysq:%s:%i: %s\n",__FILE__,__LINE__,(msg)); exit(-1);}

/**
   @brief Prints debug message specified by format fmt and a variable number of arguments
*/
#define DEBUG(fmt, ...) {printf("debug librysq:%s:%i: ", __FILE__, __LINE__); printf(fmt,__VA_ARGS__); }

/**
   @brief Creates shuffle mask 
*/
#define SHUFFLE(i,j,k,l) ((i) | (j)<<2 | (k)<<4 | (l)<<6)

namespace util {

    static inline int triangular(int n) { return (n*(n + 1))/2; }

    static inline int tetrahedral(int n) { return (n*(n + 1)*(n + 2))/6; }

    template<int n> struct binomial2 {
	static const int value = (n*n - n)/2;
    };

    template<int n> struct pow {
	template<typename T>
	static T eval(T x) { return x*pow<n-1>::eval(x); }
    };

    template<> struct pow<0> {
	template<typename T>
	static T eval(T x) { return 1; }
    };

    template<typename T>
    T pow4(T x) { return pow2(x)*pow2(x); }

    template<typename T>
    T ceiling(T m, T n) { return m/n + (m%n > 0); }

    /** 
     * Next power of two.
     * @param m 
     * @return 
     */
    template<typename T>
    T ceiling2(T n) {
	T k = 1;
	while (k < n) k *= 2;
	return k;
    }

    template<int m>
    struct round {
	template<typename T> T operator()(T n) { return n + n%m; }
    };

    struct byteset {
	typedef uint32_t type;
	typedef int8_t byte;
	typedef uint8_t ubyte;

	template<byte B0, byte B1, byte B2, byte B3>
	struct set {
	    static const type value = (type(static_cast<uint8_t>(B0)) << 0 |
				       type(static_cast<uint8_t>(B1)) << 8 |
				       type(static_cast<uint8_t>(B2)) << 16 |
				       type(static_cast<uint8_t>(B3)) << 24);
	};

	template<type B, size_t i>
	struct at { static const byte signed_value = byte((B >> i*8) & 0xFF); };

	static inline type cast(byte B) { return type(ubyte(B)); }

	static inline type value(byte B0, byte B1, byte B2, byte B3) {
	    return (cast(B0) << 0 | cast(B1) << 8 | cast(B2) << 16 | cast(B3) << 24);
	}
    };


#ifdef __CUDACC__
    static inline ushort size(const dim3 &v) { return v.x*v.y*v.z; }
#endif


    template<template <class> class Q, typename T, size_t N>
    static inline void unpack(const Q<boost::array<T,N> > &A,
			      const T* &t0, const T* &t1, const T* &t2, const T* &t3) {
	t0 = A[0].elems;
	t1 = A[1].elems;
	t2 = A[2].elems;
	t3 = A[3].elems;
    }

    template<class C, typename T>
    static inline void unpack(const C &A, T &t0, T &t1, T &t2, T &t3) {
	t0 = A[0];
	t1 = A[1];
	t2 = A[2];
	t3 = A[3];
    }

    /**
       @brief Unmasks the bit shuffling
       @param mask the bits to unmask
       @param[out] i position of i
       @param[out] j position of j
       @param[out] k position of k
       @param[out] l position of l
    */

    template <typename S, typename D> void add(int N, D *dst, const S *src) {
	for (int i = 0; i < N; ++i) {
	    dst[i] += src[i];
	}
    }

    template <typename S, typename D> void add(int N, D *dst, S a, const S *src) {
	for (int i = 0; i < N; ++i) {
	    dst[i] += a*src[i];
	}
    }

    template <typename S, typename D> void copy(int N, D *dst, const S *src) {
	for (int i = 0; i < N; ++i) {
	    dst[i] = src[i];
	}
    }

    template <int N, typename S, typename D> void copy(const S *source, D *dest) {
	for (int i = 0; i < N; ++i) {
	    dest[i] = source[i];
	}
    }

    template <typename T> inline
    void scale(int n, T a, T *A) {
	for(int i = 0; i < n; ++i) {
	    A[i] *= a;
	}
    }

    template <typename T> inline
    void scale(int n, T *dest, T s, T const *src) {
	for(int i = 0; i < n; ++i) {
	    dest[i] = s*src[i];
	}
    }

}

#undef _rysq_host_device_

#endif // _UTIL_H
