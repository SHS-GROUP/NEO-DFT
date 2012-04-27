#ifndef __CUDA_RECURRENCE_H
#define __CUDA_RECURRENCE_H

#include "rysq.hpp"
#include "cuda/kernels/device.hpp"
#include "cuda/kernels.hpp"
#include "cuda/kernels/recurrence.hpp"

using namespace device;
using namespace rysq::cuda::kernels;

// template< ushort N >__global__
// static void recurrence(int m, int n, double A, double B,
// 		       const double *rAB, const double *rAi, const double *rBk,
// 		       const double *t2, const double *W,
// 		       double *G) {
//     G[-1] = 0;
//     double A1 = 1.0/A, B1 = 1.0/B;
//     ushort bSize = block::size();
//     ushort bRank = block::rank();
//     recurrence<N>(m, n, A, B, A1, B1, rAB, rAi, rBk, t2, W, G, bSize, bRank);
// }

// // G(ij,kl,3,N) 
// // assumes G is padded with 0 element to deal with recurrence
// void rysq::cuda::kernels::recurrence(int N, int Lij, int Lkl, double A, double B, 
// 				     const double *rAB, const double *rAi, const double *rBk,
// 				     const double *t2, const double *W,
// 				     double *Gx, double *Gy, double *Gz) {

//     const int ldN = N+N%2;

// #define Gx(a,i,k) (Gx[a + ldN*(i + (Lij +1)*(k))])
// #define Gy(a,i,k) (Gy[a + ldN*(i + (Lij +1)*(k))])
// #define Gz(a,i,k) (Gz[a + ldN*(i + (Lij +1)*(k))])
// #define tG(i,k,r,a) (tG[i + (Lij +1)*(k + (Lkl+1)*(r + 3*a))])

//     double tG[(Lij+1)*(Lkl+1)*3*N];

//     dim3 gDim = dim3(1,1,1);
//     // Lij1, xyz, roots
//     dim3 bDim = dim3(8,3,N);
//     int Ns = 9*N*sizeof(double); //recurrence coeficients

//     double A1 = 1.0/A;
//     double B1 = 1.0/B;

//     double *d0 = rysq::cuda::copyToDevice(rAB, 3);
//     double *d1 = rysq::cuda::copyToDevice(rAi, 3);
//     double *d2 = rysq::cuda::copyToDevice(rBk, 3);
//     double *d3 = rysq::cuda::copyToDevice(t2, N);
//     double *d4 = rysq::cuda::copyToDevice(W, N);

//     double *gG = rysq::cuda::allocate<double>((Lij+1)*(Lkl+1)*3*N+1);

//     ::recurrence<5><<<gDim, bDim, Ns>>>(Lij+1, Lkl+1, A, B,
// 				     d0,d1,d2,d3,d4,
// 				     gG+1);

//     // compute on device
//     rysq::cuda::copyToHost(tG, gG+1, (Lij+1)*(Lkl+1)*3*N);
//     rysq::cuda::free(gG);
//     rysq::cuda::free(d0);
//     rysq::cuda::free(d1);
//     rysq::cuda::free(d2);
//     rysq::cuda::free(d3);
//     rysq::cuda::free(d4);

//     for (int k = 0; k <= Lkl; ++k) {
// 	for (int i = 0; i <= Lij; ++i) {
// 	    for (int a = 0; a < N; ++a) {
// 		Gx(a,i,k) = tG(i,k,0,a);
// 		Gy(a,i,k) = tG(i,k,1,a);
// 		Gz(a,i,k) = tG(i,k,2,a);
// 	    }
// 	}
//     }

// #undef Gx
// #undef Gy
// #undef Gz
// #undef tG

// Start}


extern __shared__ double dmem[];

// __global__
// static void recurrence(int Lij, int Lkl, double A, double B,
// 		       double A1, double B1,
// 		       const double *rAB, const double *rAi, const double *rBk,
// 		       const double *t2, const double *W,
// 		       double *G) {

//     ushort bRank = block::rank();
//     ushort bSize = block::size();
    
//     int N = blockDim.z;

//     double *Cm = dmem;
//     double *B00 = Cm + 3*N;
//     double *B10 = B00 + N;
//     double *Cn = B10 + N;
//     double *B01 = Cn + 3*N;

//     for (int a = bRank; a < N; a += bSize) {
// 	double t = t2[a];
// 	B00[a] = 0.5*t;
// 	B10[a] = 0.5*A1*(1.0 - B*t);
// 	B01[a] = 0.5*B1*(1.0 - A*t);
// 	for (int r = 0; r < 3; ++r) {
// 	    double q = rAB[r]*t;
// 	    Cm[r+a*3] = rAi[r] - B*q;
// 	    Cn[r+a*3] = rBk[r] + A*q;
// 	}
//     }

//     __syncthreads();

//     int m = Lij+1;
//     int n = Lkl+1;
//     int x = threadIdx.x;
//     int z = threadIdx.z;
//     int yz = threadIdx.y + threadIdx.z*3;
//     double *myG = G + yz*m*n;

//     if (bRank == 0)myG[-1] = 0;

//     if (x == 0) {
// 	myG[0] = (threadIdx.y == 2) ? W[z] : 1.0;
// 	myG[1] = Cm[yz]*myG[0];
// 	for(int i = 2; i < m; ++i) {
// 	    myG[i] = (i-1)*B10[z]*myG[(i-2)] + Cm[yz]*myG[(i-1)];
// 	}
//     }

//     __syncthreads();

//     if(n == 1) return;

//     if (x < m) {
// 	myG[x+m] = x*B00[z]*myG[x-1] + Cn[yz]*myG[x];
// 	for(int k = 2; k < n; ++k) {
// 	    myG += m; // next column
// 	    myG[x+m] = (k-1)*B01[z]*myG[x-m] + x*B00[z]*myG[x-1]+Cn[yz]*myG[x];
// 	}
//     }

// }


#endif // __CUDA_RECURRENCE_H
