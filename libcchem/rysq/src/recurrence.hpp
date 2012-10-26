#ifndef _RYSQ_KERNEL_RECURRENCE_HPP_
#define _RYSQ_KERNEL_RECURRENCE_HPP_


namespace rysq {

    template<int M, int N, typename T>
    static inline void recurrence(int m, T A, T B,
				  const double *__restrict rAi,  const double *__restrict rBk, 
				  const double *__restrict rAB, 
				  const double *__restrict t2, const double *__restrict W,
				  T *__restrict B0, T *__restrict B1,
				  T *__restrict Cx, T *__restrict Cy, T *__restrict Cz,
				  T *__restrict Gx, T *__restrict Gy, T *__restrict Gz) {
    
	const int ldN = N + N%2;
	T A2 = 1.0/(2*A);

	T *Gx0 = Gx;
	T *Gy0 = Gy;
	T *Gz0 = Gz;
    
	for(int a = 0; a < N; ++a) {
	    Gx0[a] = 1.0;
	    Gy0[a] = 1.0;
	    Gz0[a] = W[a];

	    T Bt2 = B*t2[a];
	
	    T cx = (rAi[0] - rAB[0]*Bt2);
	    T cy = (rAi[1] - rAB[1]*Bt2);
	    T cz = (rAi[2] - rAB[2]*Bt2);
	    Gx0[a+ldN] = cx;
	    Gy0[a+ldN] = cy;
	    Gz0[a+ldN] = cz * Gz0[a];
	    T B1 = (1.0 - Bt2)*A2;

	    for (int i = 1; i < M; ++i) {
		T di = i;
		Gx0[a+(i+1)*ldN] = di*B1*Gx0[a+(i-1)*ldN] + cx*Gx0[a+i*ldN];
		Gy0[a+(i+1)*ldN] = di*B1*Gy0[a+(i-1)*ldN] + cy*Gy0[a+i*ldN];
		Gz0[a+(i+1)*ldN] = di*B1*Gz0[a+(i-1)*ldN] + cz*Gz0[a+i*ldN];
	    }
	}

	if (m == 0) return ;

	T *Gx1 = Gx + ldN*(M+1);
	T *Gy1 = Gy + ldN*(M+1);
	T *Gz1 = Gz + ldN*(M+1);
		   
	T B2 = 1.0/(2*B);
	
	for(int a = 0; a < N; ++a) {

	    // coefficients
	
	    Cx[a] = (rBk[0] + A*rAB[0]*t2[a]);
	    Cy[a] = (rBk[1] + A*rAB[1]*t2[a]);
	    Cz[a] = (rBk[2] + A*rAB[2]*t2[a]);
	    B0[a] = 0.5*t2[a];
	    B1[a] = (1.0 - A*t2[a])*B2;
	
	    Gx1[a] = Cx[a]*Gx0[a];
	    Gy1[a] = Cy[a]*Gy0[a];
	    Gz1[a] = Cz[a]*Gz0[a];
	
	    for (int i = 1; i <= M; ++i) {
		T di = i;
		Gx1[a+i*ldN] = di*B0[a]*Gx0[a+(i-1)*ldN] + Cx[a]*Gx0[a+i*ldN];
		Gy1[a+i*ldN] = di*B0[a]*Gy0[a+(i-1)*ldN] + Cy[a]*Gy0[a+i*ldN];
		Gz1[a+i*ldN] = di*B0[a]*Gz0[a+(i-1)*ldN] + Cz[a]*Gz0[a+i*ldN];
	    }
	}
    }

    template<int M, int N, typename T>
    static inline void recurrence(int m, T *B0, T *__restrict B1,
				  T *__restrict Cx, T *__restrict Cy, T *__restrict Cz,
				  T *__restrict Gxm1, T *__restrict Gxm2, T *__restrict Gxm,
				  T *__restrict Gym1, T *__restrict Gym2, T *__restrict Gym,
				  T *__restrict Gzm1, T *__restrict Gzm2, T *__restrict Gzm) {
 
	static const int ldN = N + N%2;

	T m1 = m - 1;
    
	for(int a = 0; a < N; ++a) {
	    Gxm[a] = m1*B1[a]*Gxm2[a] + Cx[a]*Gxm1[a];
	    Gym[a] = m1*B1[a]*Gym2[a] + Cy[a]*Gym1[a];
	    Gzm[a] = m1*B1[a]*Gzm2[a] + Cz[a]*Gzm1[a];
	
	    for (int i = 1; i <= M; ++i) {
		T di = i;
		Gxm[a+i*ldN] = m1*B1[a]*Gxm2[a+i*ldN] + di*B0[a]*Gxm1[a+(i-1)*ldN] + Cx[a]*Gxm1[a+i*ldN];
		Gym[a+i*ldN] = m1*B1[a]*Gym2[a+i*ldN] + di*B0[a]*Gym1[a+(i-1)*ldN] + Cy[a]*Gym1[a+i*ldN];
		Gzm[a+i*ldN] = m1*B1[a]*Gzm2[a+i*ldN] + di*B0[a]*Gzm1[a+(i-1)*ldN] + Cz[a]*Gzm1[a+i*ldN];
	    }
	}
    }

    template<int M, int N, typename T>
    static void recurrence(int m, double A, double B,
			   const double *__restrict rAB,
			   const double *__restrict rAi, const double *__restrict rBk, 
			   const double *__restrict t2, const double *__restrict W,
			   T *__restrict Gx, T *__restrict Gy, T *__restrict Gz) {

	static const int N16 = N + N%2;
#define Gx(j) (&Gx[(j)*(N16*(M+1))])
#define Gy(j) (&Gy[(j)*(N16*(M+1))])
#define Gz(j) (&Gz[(j)*(N16*(M+1))])

	T B0[N16] __attribute__ ((aligned(16)));
	T B1[N16] __attribute__ ((aligned(16)));
	T Cx[N16] __attribute__ ((aligned(16)));
	T Cy[N16] __attribute__ ((aligned(16)));
	T Cz[N16] __attribute__ ((aligned(16)));
    
	recurrence<M,N>(m, T(A), T(B), rAi, rBk, rAB, t2, W, 
			B0, B1, Cx, Cy, Cz, Gx, Gy, Gz);
	
	for(int j = 2; j <= m; ++j) {
	    recurrence<M,N>(j, B0, B1, Cx, Cy, Cz, 
			    Gx(j-1), Gx(j-2), Gx(j), 
			    Gy(j-1), Gy(j-2), Gy(j),
			    Gz(j-1), Gz(j-2), Gz(j));
	}

#undef Gx
#undef Gy
#undef Gz
    
    }

}


#endif /* _RYSQ_KERNEL_RECURRENCE_HPP_ */
