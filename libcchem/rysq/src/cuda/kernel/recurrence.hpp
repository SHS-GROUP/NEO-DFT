#ifndef RYSQ_CUDA_KERNEL_RECURRENCE_HPP
#define RYSQ_CUDA_KERNEL_RECURRENCE_HPP


template<typename T>
__device__
void recurrence1(int N, int m, int n,
		 const double &A_, const double &B,
		 const double *ri, const double *rk,
		 const double *rA, const double *rB,
		 const double *t2, const double *W,
		 T* G, const int3 &thread) {
    if (thread.z < N) {
	const int &y = thread.y;
	const int &z = thread.z;
	T *g = G + (y + z*3)*m*n;
	T t = t2[z];
    	T b = 0.5*(1.0 - B*t)*A_;///(2*A);
    	T Cm = (rA[y] - ri[y]) - t*B*(rA[y] - rB[y]);
	T g0 = (y < 2) ? 1 : W[z];
	T g1 = Cm*g0;
    	g[0] = g0;
    	g[1] = g1;
    	for (int i = 2; i < m; ++i) {
	    T gi = (i-1)*b*g0 + Cm*g1;
	    g0 = g1;
	    g1 = gi;
    	    g[i] = gi;//(i-1)*b*g[i-2] + Cm*g[i-1];
	    //g[i] = (i-1)*b*g[i-2] + Cm*g[i-1];
    	}
    }
}

/**
   n must not be constant expression due to compiler bug
 */
template<typename T>
__device__
void recurrence2(int N, int m, int n,
		 const double &A, const double &B_,
		 const double *ri, const double *rk,
		 const double *rA, const double *rB,
		 const double *t2, const double *W,
		 T* G, const int3 &thread) {

    if (n == 1) return;

    if (thread.y < 3 && thread.z < N && thread.x < m) {
	const T &t = t2[thread.z];
	T Cn = (rB[thread.y]-rk[thread.y]) + A*(rA[thread.y] - rB[thread.y])*t;
	T B00 = thread.x*0.5*t;//B00[thread.z];

	T *g = G + thread.x + (thread.y + thread.z*3)*m*n;
	const T *g0 = g-m;

	g[m] = ((thread.x) ? B00*g[-1] : 0) + Cn*g[0]; // avoid G[-1]

	//return;
	if (n > 2) {
	    T gk0 = g[0*m];
	    T gk1 = g[m];
	    T B01 = 0.5*(1.0 - A*t)*B_;///(2*B);//B01/*[thread.z]*/;
	    {
		int k = 2;
		T gk2 = (k-1)*B01*gk0 /*+ B00*g[(k-1)*m-1]*/ + Cn*gk1;
		gk2 += B00*g0[k*m-1];
		gk0 = gk1;
		gk1 = gk2;
		g[k*m] = gk2;
	    }
#pragma unroll 1
	    for (int k = 3; k < n; ++k) {
		T gk2 = (k-1)*B01*gk0 /*+ B00*g[(k-1)*m-1]*/ + Cn*gk1;
		gk2 += B00*g0[k*m-1];
		gk0 = gk1;
		gk1 = gk2;
		g[k*m] = gk2;
	    }
	}

    }

}

template<int N, typename T>
__device__
void recurrence(const ushort &m, const ushort &n,
		const double &A, const double &B,
		const double &rA, const double &rB,
		const double &rAi, const double &rBk,
	        const double &t2, const double &W,
		T* G,
		const dim3 &thread) { 

    T *myG = G + (thread.y + thread.z*3)*m*n;

    if (thread.x == 0 && thread.y < 3 && thread.z < N) {
    	const T &t = t2/*[thread.z]*/;
    	T b = (1.0 - B*t)/(2*A);
    	T Cm = rAi/*[thread.y]*/ - t*B*(rA/*[thread.y]*/ - rB/*[thread.y]*/);
    	myG[0] = (thread.y < 2) ? 1 : W/*[thread.z]*/;
    	myG[1] = Cm*myG[0];
    	for (ushort i = 2; i < m; ++i) {
    	    myG[i] = (i-1)*b*myG[i-2] + Cm*myG[i-1];
    	}
    }

    if (n == 1) return;

    if (thread.x < m && thread.y < 3 && thread.z < N) {
	const T t = t2/*[thread.z]*/;
	T Cn = rBk/*[thread.y]*/ + A*(rA/*[thread.y]*/ - rB/*[thread.y]*/)*t;
	T B00 = thread.x*0.5*t;//B00/*[thread.z]*/;

	T *g = myG + thread.x;
	const T *g0 = g-m;

	g[m] = ((thread.x) ? B00*g[-1] : 0) + Cn*g[0]; // avoid G[-1]

	//return;
	if (n == 2) return;

	T gk0 = g[0*m];
	T gk1 = g[m];
	T B01 = (1.0 - A*t)/(2*B);//B01/*[thread.z]*/;

#pragma unroll 1
	for (ushort k = 2; k < n; ++k) {
	    T gk2 = (k-1)*B01*gk0 /*+ B00*g[(k-1)*m-1]*/ + Cn*gk1;
	    gk2 += B00*g0[k*m-1];
	    gk0 = gk1;
	    gk1 = gk2;
	    g[k*m] = gk2;
	}

    }

}


#endif /* RYSQ_CUDA_KERNEL_RECURRENCE_HPP */
