#ifndef RYSQ_ROOTS_EVALUATE_HPP
#define RYSQ_ROOTS_EVALUATE_HPP

#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/repetition/enum_binary_params.hpp>

// namespace rysq {
// namespace detail {


    template<typename U>
    BOOST_GPU_ENABLED
    inline bool master(const U &thread){ return thread == 0; }

    BOOST_GPU_ENABLED
    inline bool master(const serial_tag&) { return true; }

    template<typename U>
    BOOST_GPU_ENABLED
    inline const U& index(const U &thread){ return thread; }

    BOOST_GPU_ENABLED
    inline size_t index(const serial_tag&) { return 0; }

    template<int L, typename T>
    BOOST_GPU_ENABLED
    void evaluate1(T X, T &t2, const T (&C)[L]) {
	T q = X;
	t2 = C[0] + C[1]*q;
	for (int i = 2; i < L; ++i) {
	    q *= X;
	    t2 += C[i]*q;
	}
    }


    template<int U>
    struct unroll {

	template<int I, class T = void>
	struct int_ {
	    typedef int enabled;
	    BOOST_GPU_ENABLED operator int() const { return I; }
	};
	template<class T>
	struct int_<U,T> {};

#define UNROLL_APPLY(N)							\
	template<class F, BOOST_PP_ENUM_PARAMS(N, class A)>		\
	BOOST_GPU_ENABLED							\
	static void apply(const F &f,					\
			  BOOST_PP_ENUM_BINARY_PARAMS(N, A, &a)) {	\
	    apply(f, BOOST_PP_ENUM_PARAMS(N, a), int_<0>());		\
	}								\
	template<class F, BOOST_PP_ENUM_PARAMS(N, class A), int I>	\
	BOOST_GPU_ENABLED							\
	static void apply(const F &f,					\
			  BOOST_PP_ENUM_BINARY_PARAMS(N, A, &a),	\
			  int_<I>,					\
			  typename int_<I>::enabled = 0) {		\
	    f(BOOST_PP_ENUM_PARAMS(N, a), int_<I>());			\
	    apply(f, BOOST_PP_ENUM_PARAMS(N, a), int_<I+1>());		\
	}								\
	template<class F, BOOST_PP_ENUM_PARAMS(N, class A)>		\
	BOOST_GPU_ENABLED							\
	static void apply(const F &f,					\
			  BOOST_PP_ENUM_BINARY_PARAMS(N, A, &a),	\
			  int_<U>) {}
	UNROLL_APPLY(3)
	UNROLL_APPLY(4)
	UNROLL_APPLY(5)

#undef UNROLL_APPLY
    };


    namespace polynomial {

	namespace functor {

	    struct evaluate {
		template<int N, int L, class A, typename T, typename U>
		BOOST_GPU_ENABLED
		void operator()(A &R, const T X, const T (&C)[N][L],
				const U &thread) const {
		    T q = X;
		    R[thread] = C[thread][0] + C[thread][1]*q;
		    for (int i = 2; i < L; ++i) {
			q *= X;
			R[thread] += C[thread][i]*q;
		    }
		}
	    };

	    struct evaluate2 {
		template<int N, int L, int M, class A, typename T, typename U>
		BOOST_GPU_ENABLED
		void operator()(const T &X, A &R, A &W,
				const T (&C)[N][L], const T (&D)[N][M],
				const U &thread) const {
		    T q = X;
		    R[thread] = C[thread][0] + C[thread][1]*q;
		    W[thread] = D[thread][0] + D[thread][1]*q;
		    for (int i = 2; i < L; ++i) {
			q *= X;
			R[thread] += C[thread][i]*q;
			W[thread] += D[thread][i]*q;
		    }
		    for (int i = L; i < M; ++i) {
			q *= X;
			W[thread] += D[thread][i]*q;
		    }
		}
	    };

	    struct inverse {
		template<int N, int L, class A, typename T, typename U>
		BOOST_GPU_ENABLED
		void operator()(A &R, const T (&C)[N][L], T x,
				const U &thread) const {
		    T p = 0.0;
		    for (int i = L - 1; i > 0; --i) {
			p = (p + C[thread][i])/x;
		    }
		    R[thread] = (C[thread][0] + p);
		}
		template<int L, typename T>
		BOOST_GPU_ENABLED
		void operator()(T &R, const T (&C)[L], T x) const {
		    R = 0.0;
		    for (int i = L - 1; i > 0; --i) {
			R = (R + C[i])/x;
		    }
		    R += C[0];
		    // return (C[0] + p);
		}
	    };

	    struct inverse0 {
		template<int N, int L, class A, typename T, typename U>
		BOOST_GPU_ENABLED
		void operator()(A &R, const T (&C)[N][L], T x,
				const U &thread) const {
		    T p = 0.0;
		    for (int i = L - 1; i >= 0; --i) {
			p = (p + C[thread][i])/x;
		    }
		    R[thread] += p;
		}
	    };

	}

	template<int N, int L, class A, typename T, typename U>
	BOOST_GPU_ENABLED
	void evaluate(A &R, const T X, const T (&C)[N][L],
		      const U &thread) {
	    functor::evaluate()(R, X, C, thread);
	}

	template<int N, int L, int M, class A, typename T, typename U>
	BOOST_GPU_ENABLED
	void evaluate2(T X, A &R, A &W,
		       const T (&C)[N][L], const T (&D)[N][M],
		       const U &thread) {
	    functor::evaluate2()(X, R, W, C, D, thread);
	}


	template<int N, typename T>
	BOOST_GPU_ENABLED
	T inverse(const T (&C)[N], T x) {
	    T R;
	    functor::inverse()(R, C, x);
	    return R;
	    
	}

	template<int N, int L, class A, typename T, typename U>
	BOOST_GPU_ENABLED
	void inverse(A &R, const T (&C)[N][L], T x,
		     const U &thread) {
	    functor::inverse()(R, C, x, thread);
	}

	template<int N, int L, class A, typename T, typename U>
	BOOST_GPU_ENABLED
	void inverse0(A &R, const T (&C)[N][L], T x,
		      const U &thread) {
	    functor::inverse0()(R, C, x, thread);
	}

	template<int N, int L, typename T>
	BOOST_GPU_ENABLED
	void evaluate(T (&R)[N], T x, const T (&C)[N][L], 
		      const serial_tag& = serial_tag()) {
	    unroll<N>::apply(functor::evaluate(), R, x, C);
	}

	template<int N, int L, int M, typename T>
	BOOST_GPU_ENABLED
	void evaluate2(T X, T (&R)[N], T (&W)[N],
		       const T (&C)[N][L], const T (&D)[N][M],
		       const serial_tag& = serial_tag()) {
	    unroll<N>::apply(functor::evaluate2(), X, R, W, C, D);
	}


	template<int N, int L, typename T>
	BOOST_GPU_ENABLED
	void inverse(T (&R)[N], const T (&C)[N][L], T x,
		     const serial_tag& = serial_tag()) {
	    unroll<N>::apply(functor::inverse(), R, C, x);
	}

	template<int N, int L, typename T>
	BOOST_GPU_ENABLED
	void inverse0(T (&R)[N], const T (&C)[N][L], T x,
		      const serial_tag& = serial_tag()) {
	    unroll<N>::apply(functor::inverse0(), R, C, x);
	}



    }

    template<int N, typename T>
    BOOST_GPU_ENABLED
    T evaluate_polynomial(const T (&C)[N], T x) {
	// evaluate1(x, t2, C, a);
	T p = 0;
// #pragma unroll
	for (int i = N-1; i > 0; --i) {
	    p = (p + C[i])*x;
	}
	return (p + C[0]);
    }

    template<int L, int M, typename T>
    BOOST_GPU_ENABLED
    void evaluate1(T X, T &R, T &W,
		   const T (&C)[L], const T (&D)[M]) {
	T q = X;
	R = C[0] + C[1]*q;
	W = D[0] + D[1]*q;
	for (int i = 2; i < L; ++i) {
	    q *= X;
	    R += C[i]*q;
	    W += D[i]*q;
	}
	for (int i = L; i < M; ++i) {
	    q *= X;
	    W += D[i]*q;
	}
    }

    template<int N, class A, typename T, typename U>
    BOOST_GPU_ENABLED
    void scale_add(const T e, A &R, A &W,
		   const T X, const T (&RN)[N],
		   const T w, const T (&WN)[N],
		   const U &thread) {
	R[thread] = e*R[thread] + RN[thread]/(X - RN[thread]);
	W[thread] = e*W[thread] + w*WN[thread];
    }

    template<int N, typename T>
    BOOST_GPU_ENABLED
    void scale_add(const T e, T (&R)[N], T (&W)[N],
		   const T X, const T (&RN)[N],
		   const T w, const T (&WN)[N],
		   const serial_tag& = serial_tag()) {
	for (int a = 0; a < N; ++a) scale_add(e, R, W, X, RN, w, WN, a);
    }

    template<int N, typename T>
    BOOST_GPU_ENABLED
    void scale_add(const T e1, T (&R)[N], const T e2, T (&W)[N],
		   const T X, const T (&RN)[N],
		   const T w, const T (&WN)[N],
		   const serial_tag& = serial_tag()) {
	for (int a = 0; a < N; ++a) {
	    R[a] = e1*R[a] + RN[a]/(X - RN[a]);
	    W[a] = e2*W[a] + w*WN[a];
	}
    }

    template<int N, class A, typename T, typename U>
    BOOST_GPU_ENABLED
    void scale_add(const T e, A &R, const T X, const T (&RN)[N],
		   const U &thread) {
	R[thread] = e*R[thread] + RN[thread]/(X - RN[thread]);
    }

    template<int N, typename T>
    BOOST_GPU_ENABLED
    void scale_add(const T e, T (&R)[N], const T X, const T (&RN)[N],
		   const serial_tag& = serial_tag()) {
	for (int a = 0; a < N; ++a) scale_add(e, R, X, RN, a);
    }

    template<typename T>
    BOOST_GPU_ENABLED
    void scale_add(const T e, T &R, T &W,
		   const T X, const T RN, const T w, const T WN) {
	R = e*R + RN/(X - RN);
	W = e*W + w*WN;
    }

    template<typename T>
    BOOST_GPU_ENABLED
    void scale_add(const T e1, T &R, const T e2, T &W,
		   const T X, const T (&RN), const T w, const T (&WN)) {
	R = e1*R + RN/(X - RN);
	W = e2*W + w*WN;
    }

    template<int N, typename T>
    BOOST_GPU_ENABLED
    T sum(const T *W) {
	T q = W[N-1];
	for (int i = N-2; i >= 0; --i) {
	    q += W[i];
	}
	return q;
    }

    template<int N, class A, typename T, typename U>
    BOOST_GPU_ENABLED
    void evaluate_inf(A &R, T r, const T (&CR)[N], A &W, T w, const T (&CW)[N],
		      const U &thread, bool simt = true) {
	R[thread] = CR[thread]/(r - CR[thread]);
	W[thread] = w*CW[thread];
	if (simt && master(thread)) W[0] -= sum<N-1>(W+1);
    }

    template<int N, typename T>
    BOOST_GPU_ENABLED
    void evaluate_inf(T (&R)[N], T r, const T (&CR)[N],
		      T (&W)[N], T w, const T (&CW)[N],
		      const serial_tag& = serial_tag()) {
	for (int i = 0; i < N; ++i) evaluate_inf(R, r, CR, W, w, CW, i, false);
	W[0] -= sum<N-1>(W+1);
    }

    template<typename T>
    BOOST_GPU_ENABLED
    void change_variable(T &R) {
	R = R/(1.0 + R);
    }

    template<int N, class A, typename U>
    BOOST_GPU_ENABLED
    void change_variable(A &R, const U &thread) {
	R[thread] = R[thread]/(1.0 + R[thread]);
    }

    template<int N, typename T>
    BOOST_GPU_ENABLED
    void change_variable(T (&R)[N], const serial_tag& = serial_tag()) {
	for(int i = 0; i < N; ++i) change_variable<N>(R,i);
    }


    template<typename T>
    BOOST_GPU_ENABLED
    void scale(T &p, const T &x, const T &r){
	p = x*r;
    }

    template<typename T>
    BOOST_GPU_ENABLED
    void divide(T &p, const T &x, const T &r){
	p = r/x;
    }

    template<size_t N, class A, typename T, typename U>
    BOOST_GPU_ENABLED
    void divide(A &p, const T &x, const T (&r)[N], const U &thread) {
    	p[thread] = r[thread]/x;
    }

    template<size_t N, class A, typename T, typename U>
    BOOST_GPU_ENABLED
    void scale(A &p, const T &x, const T (&r)[N], const U &thread) {
    	p[thread] = r[thread]*x;
    }

    template<size_t N, typename T>
    BOOST_GPU_ENABLED
    void divide(T (&p)[N], const T &x, const T (&r)[N],
		const serial_tag &tag = serial_tag()) {
    	for (size_t i = 0; i < N; ++i) divide(p, x, r, i);
    }

    template<size_t N, class A, typename T, typename U>
    BOOST_GPU_ENABLED
    inline void scale(T s, A &R, const U &thread) {
	R[thread] *= s;
    }

    template<size_t N, typename T>
    BOOST_GPU_ENABLED
    void scale(T (&p)[N], T x, const T (&r)[N],
	       const serial_tag &tag = serial_tag()) {
    	for (size_t i = 0; i < N; ++i) scale(p, x, r, i);
    }

    template<size_t N, typename T>
    BOOST_GPU_ENABLED
    void scale(T s, T (&R)[N], const serial_tag& = serial_tag()) {
	for (size_t i = 0; i < N; ++i) scale<N>(s, R, i);
    }

    template<int M, int N, class A, typename T, typename U>
    BOOST_GPU_ENABLED
    void add_evaluate_polynomial1(A &R, const T (&C)[M][N], T x,
				  const U &thread) {
	T p = 0.0;
	for (int i = N - 1; i >= 0; --i) {
	    p = (p + C[thread][i])*x;
	}
	R[thread] += p;
    }

    template<int M, int N, typename T>
    BOOST_GPU_ENABLED
    void add_evaluate_polynomial1(T (&R)[M], const T (&C)[M][N], T x,
				  const serial_tag&) {
	for (int i = 0; i < M; ++i)
	    add_evaluate_polynomial1(R, C, x, i);
    }


    template<int N, class A, typename T, typename U>
    BOOST_GPU_ENABLED
    void add_evaluate_inf(A &R, const T (&CR)[N], T r, const U &thread) {
	R[thread] += CR[thread]/(r - CR[thread]);
    }

    template<int N, typename T>
    BOOST_GPU_ENABLED
    void add_evaluate_inf(T (&R)[N], const T (&CR)[N], T r,
			  const serial_tag& = serial_tag()) {
	// #pragma unroll
	for (int i = 0; i < N; ++i)
	    add_evaluate_inf(R, CR, r, i);
    }


// } // namespace detail
// } // namespace rysq

#endif /* RYSQ_ROOTS_EVALUATE_HPP */

