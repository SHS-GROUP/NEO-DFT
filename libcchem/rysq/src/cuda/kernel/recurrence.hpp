#ifndef _RYSQ_CUDA_KERNEL_RECURRENCE_HPP_
#define _RYSQ_CUDA_KERNEL_RECURRENCE_HPP_

template<int N, typename T> __device__
void recurrence(const ushort &m, const ushort &n,
		const double &A, const double &B,
		const double &A1, const double &B1,
		const double *rAB, const double *rAi, const double *rBk,
	        const double *t2, const double *W, T *G,
		const ushort &bSize, const ushort &bRank) {

//     using namespace device;
//     ushort bRank = block::rank();
//     ushort bSize = block::size();

    if (bRank == 0) {
	ushort mn = m*n;
	for (ushort i = 0, j = 0; i < N; ++i) {
	    G[j] = 1.0;
	    j += mn;
	    G[j] = 1.0;
	    j += mn;
	    G[j] = W[i];
	    j += mn;
	}
    }

    __syncthreads();

    __shared__ T Cm[N*3], Cn[N*3];
    __shared__ T B00[N], B10[N], B01[N];

    //for (ushort a = bRank; a < N; a += bSize) {
    if (bRank < N) {
	const ushort &a = bRank;
	T t = t2[a];
	B00[a] = 0.5*t;
	B10[a] = 0.5*A1*(1.0 - B*t);
	B01[a] = 0.5*B1*(1.0 - A*t);
	for (ushort r = 0; r < 3; ++r) {
	    T q = rAB[r]*t;
	    Cm[r+a*3] = rAi[r] - B*q;
	    Cn[r+a*3] = rBk[r] + A*q;
	}
    }

    __syncthreads();

    short x = threadIdx.x;
    const ushort &y = threadIdx.y;
    const ushort &z = threadIdx.z;
    ushort yz = threadIdx.y + threadIdx.z*3;
    T *myG = G + yz*m*n;

    if (x == 0 and z < N) {
	T B = B10[z];
	T C = Cm[yz];
	myG[1] = C*myG[0];

	for(ushort i = 2; i < m; ++i) {
	    myG[i] = B*myG[i-2] + C*myG[i-1];
	    B += B10[z];
	}
    }

    __syncthreads();

    if(n == 1) return;

    if (x < m and z < N) {
	T C = Cn[yz];
	T B0 = x*B00[z];
	myG[x+m] = ((x) ? B0*myG[x-1] : 0) + C*myG[x]; // avoid G[-1]
	T q = B01[z];
	for(ushort k = 2; k < n; ++k) {
	    myG += m; // next column
	    myG[x+m] = q*myG[x-m] + B0*myG[x-1] + C*myG[x];
	    q += B01[z];
	}
    }

}


#endif /* _RYSQ_CUDA_KERNEL_RECURRENCE_HPP_ */
