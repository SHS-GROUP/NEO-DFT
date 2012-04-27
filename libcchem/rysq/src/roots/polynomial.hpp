// #ifndef _RYSQ_ROOTS_POLYNOMIAL_H_
// #define _RYSQ_ROOTS_POLYNOMIAL_H_

// namespace {

//     template<int N, int n> __device__
//     static inline void polynomial(double *p, const double x, const double (&c)[N][n]) {
// 	for (int a = 0; a < N; ++a) {
// 	    p[a] = c[a][0];
// 	    double q = 1.0;
// 	    for (int i = 1; i < n; ++i) {
// 		q *= x;
// 		p[a] += q*c[a][i];
// 	    }
// 	}
//     }

//     template<int N, int n> __device__
//     static inline void polynomial(double *p, const double x, const double (&c)[N][n],
// 		    unsigned short a) {
// 	p[a] = c[a][0];
// 	double q = 1.0;
// 	for (int i = 1; i < n; ++i) {
// 	    q *= x;
// 	    p[a] += q*c[a][i];
// 	}
//     }

//     template<int N, int n> __device__
//     static inline void add_negative(double *p, const double x, const double (&c)[N][n]) {
// 	for (int a = 0; a < N; ++a) {
//  	    double q = 0.0;
// // 	    for (int i = n-1; i >= 0; --i) {
// // 		q = (q + c[a][i])/x;
// // 	    }
// // 	    p[a] += q;
// 	    double y = x;
//  	    for (int i = 0; i < n-1; ++i) {
// 		q = (q + c[a][i])*x;
// 		y *= x;
// 	    }
// 	    q += c[a][n-1];
// 	    p[a] += q/y;
// 	}
//     }

//     template<int N, int n> __device__
//     static inline void add_negative(double *p, const double x, const double (&c)[N][n],
// 		      unsigned short a) {
// 	double q = 0.0;
// // 	for (int i = n-1; i >= 0; --i) {
// // 	    q = (q + c[a][i])/x;
// // 	}
// // 	p[a] += q;
// 	    double y = x;
//  	    for (int i = 0; i < n-1; ++i) {
// 		q = (q + c[a][i])*x;
// 		y *= x;
// 	    }
// 	    q += c[a][n-1];
// 	    p[a] += q/y;
//     }

//     template<int N> __device__
//     static inline void scale_add(const double e, double *p, const double x, const double r[N]){
// 	for (int a = 0; a < N; ++a) {
// 	    p[a] = e*p[a] + x*r[a];
// 	}
//     }

//     template<int N> __device__
//     static inline void scale_add(const double e, double *p, const double x, const double r[N],
// 		   unsigned short a){
// 	p[a] = e*p[a] + x*r[a];
//     }

//     template<int N> __device__
//     static inline void scale_add2(const double e, double *p, const double x, const double (&r)[N]){
// 	for (int a = 0; a < N; ++a) {
// 	    p[a] = e*p[a] + r[a]/(x - r[a]);
// 	}
//     }

//     template<int N> __device__
//     static inline void scale_add2(const double e, double *p, const double x, const double (&r)[N],
// 		    unsigned short a){
// 	p[a] = e*p[a] + r[a]/(x - r[a]);
//     }



//     template<int N> __device__
//     static inline void divide(double *p, const double x, const double r[N], unsigned short a){
// 	    p[a] = r[a]/x;
//     }



//     template<int N> __device__
//     static inline void scale(double *p, const double x, const double r[N], unsigned short a) {
// 	p[a] = x*r[a];
//     }



// //     template<int N, int n> __device__
// //     double weight0(const double X, const double (&C)[n], const double e,
// // 		   const double *W) {
// // 	double p = 0.0;
// // // 	double q = X;
// // // 	for (int i = 1; i < n-1; ++i) {
// // // 	    p = (p + C[i])*X;
// // // 	    q *= X;
// // // 	}
// // // 	p = (p + C[n-1])/q;
// // 	for (int i = n - 1; i > 0; --i) {
// // 	    p = (p + C[i])/X;
// // 	}
// // 	return e*(p + C[0]) + W[0] - sum<N-1>(&W[1]);
// //     }



//     template<int N> __device__
//     static inline void change_variable(double *t2, unsigned short a) {
// 	t2[a] = t2[a]/(1.0 + t2[a]);
//     }

// }

// #endif /* _RYSQ_ROOTS_POLYNOMIAL_H_ */
