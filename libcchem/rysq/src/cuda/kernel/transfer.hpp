#ifndef RYSQ_CUDA_KERNEL_TRANSFER_HPP
#define RYSQ_CUDA_KERNEL_TRANSFER_HPP

#include <cuda.h>
#include <cuda_runtime.h>

#include <boost/utility/binary.hpp>
#include <boost/utility/enable_if.hpp>

namespace rysq {
namespace cuda {
namespace kernel {

template<typename T, typename U>
__device__
void transfer1(int N, int mi, int mj, int n,
	       const double *dr,
	       T* __restrict__ G, T* __restrict__ I,
	       const U &thread) {

    // using namespace device;

    const int &x = thread.x;
    const int &y = thread.y;
    const int &z = thread.z;

    if (z >= N) return;
    if (y >= 3) return;

    const int mij = mi*mj;
    const int m = mi + mj - 1;

    const T dq = dr[y];
    const int yz = (y + z*3)*n;

    for (int k = 0; k < n; ++k) {
    	T *g = G + (k + yz)*(m);
    	T *I_ = I + (k + yz)*(mij);
    	if (x < mi) {
    	    I_[x] = g[x];
    	}
#pragma unroll 1
    	for (int j = 1; j < mj; ++j) {
	    I_ += mi;
	    double q = 0;
    	    if (x < m-j) {
    		q = dq*g[x] + g[x+1];
		g[x] = q;
    	    }
    	    if (x < mi) {
    		I_[x] = q;//g[x];
    	    }
    	}
    }
    return;

    //for(int kl = x; kl < n; kl += blockDim.x) {
    if (x < n) {
	const int &kl = x;
	T *myG = G + (kl + yz)*(m);
	T *myI = I + (kl + yz)*(mij);
	// printf("%p %u %u %u\n", myI, z, y, n);
 	for (int i = 0; i < mi; ++i) {
	    myI[i] = myG[i];
	}
	for (int j = 1; j < mj; ++j) {
	    T q0 = myG[0];
	    myI += mi;
	    for (int i = 0; i < mi; ++i) {
		//T q0 = myG[i];
		T q1 = myG[i+1];
		q0 = dq*q0 + q1;
		myI[i] = q0;
		myG[i] = q0;
		q0 = q1;
	    }
	    for (int i = mi; i < (m - j); ++i) {
		//myG[i] = dq*myG[i] + myG[i+1];
		T q1 = myG[i+1];
		q0 = dq*q0 + q1;
		myG[i] = q0;
		q0 = q1;
	    }
	}
    }

}

template<typename T, typename U, typename V>
__device__
void transfer2(int N, int mij, int nk, int nl,
	       const double *dr,
	       T* __restrict__ G, T* __restrict__ I,
	       const U &thread, const V &dim) {

    const int &x = thread.x;
    const int &y = thread.y;
    const int &z = thread.z;

    if (z >= N) return;
    if (y >= 3) return;
    if (x >= dim.x) return;

    const int nkl = nk*nl;
    const int n = nk + nl - 1;

    const int yz = (y + z*3)*mij;
    const T dq = dr[y];

    for (int ij = x; ij < mij; ij += dim.x) {
	T *myG = G + ij + yz*n;
	T *(myI) = I + ij + yz*nkl;

	for (int k = 0; k < nk; ++k) {
	    *(myI) = *(myG);
	    // printf("(%i,%i,%i) k=%i p=%p dq=%e v=%e\n",
	    //        ij, y, z, k, myI, dq, *(myI) );
	    myI += mij;
	    myG += mij;
	}
	myG -= nk*mij;
	for (int l = 1; l < nl; ++l) {
	    T q0 = *(myG);
	    T *Gl = myG;
	    for (int k = 0; k < nk; ++k) {
		T q1 = *(Gl + mij);
		q0 = dq*q0 + q1;
		*(myI) = q0;
		//printf("(%i,%i,%i) k=%i l=%i p=%p\n", ij, y, z, k, l, myI);
		//printf("y=%i p=%p dq=%e, q=%e\n", y, myI, dq, q0);
		*(Gl) = q0;
		q0 = q1;
		Gl += mij;
		myI += mij;
	    }
	    for (int k = nk; k < (n - l); ++k) {
		T q1 = *(Gl + mij);
		*(Gl) = dq*q0 + q1;
		q0 = q1;
		Gl += mij;
	    }
	}
    }
}


}
}
}

#endif /* RYSQ_CUDA_KERNEL_TRANSFER_HPP */
