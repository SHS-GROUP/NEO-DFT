#ifndef RYSQ_KERNEL_TRANSFER_HPP_
#define RYSQ_KERNEL_TRANSFER_HPP_

#include <assert.h>
#include <math.h>

namespace rysq {
namespace kernel {

template<int A, int B, int N, template<typename, size_t> class align>
struct transfer {

    template<typename T, typename U>
    static void apply(int C, int D, U dij, U dkl, U *__restrict G,
		      T *__restrict I, U *__restrict mem16) {
	if (B + D == 0) return;
	else if (D == 0) apply(C, dij, G, I);
	else if (B == 0) apply<(A+1)>(C, D, dkl, G, I);
	else {
	    apply(C+D, dij, G, mem16);
	    apply<(A+1)*(B+1)>(C, D, dkl, mem16, I);
	}
    }

    template<typename T, typename U>
    static void apply(int CD, U dq, U *__restrict G, T *__restrict I) {

	static const int NT = N + align<T,N>::value;
	static const int NU = N + align<U,N>::value;

#define G(n,ab) G[(n) + NU*((ab))]
#define I(n,a,b) I[(n) + NT*((a)+(A+1)*((b)))]

	for(int kl = 0; kl <= CD; ++kl) {
	
#ifdef __INTEL_COMPILER
#pragma ivdep
#pragma vector aligned
#endif
	    for(int a = 0; a < N; ++a) {
		for (int i = 0; i <= A; ++i) {
		    I(a,i,0) = G(a,i);
		}

		for (int j = 1; j <= B; ++j) {
		    for (int i = 0; i <= A; ++i) {
			G(a,i) = dq*G(a,i) + G(a,i+1);
			I(a,i,j) = G(a,i);
		    }
		    for (int i = A+1; i <= A+B-j; ++i) {
			G(a,i) = dq*G(a,i) + G(a,i+1);
		    }
		}
	    }
	    G += NU*(A+B+1);
	    I += NT*(A+1)*(B+1);
	}

#undef G
#undef I
    }
    
    template<int M, typename T, typename U>
    //  __restrict qualifier on G may cause incorrect code
    static void apply(int C, int D, U dq, U *G, T *__restrict I) {

	static const int NU = N + align<U,N>::value;
	static const int NT = N + align<T,N>::value;

	static const int MU = M*NU;
	static const int MT = M*NT;

#define G(n,ij) G[(n) + MU*((ij))]
#define I(n,i,j) I[(n) + MT*((i)+(C+1)*((j)))]

	for (int i = 0; i <= C; ++i) {
	    T * __restrict I_ = &I(0,i,0);
	    U *G_ = &G(0,i);
#ifdef __IMTEL_COMPILER
#pragma ivdep
#pragma vector aligned
#endif
	    for (int ij = 0; ij < M; ++ij) {
	    	for(int a = 0; a < N; ++a) {
	    	    I_[a] = G_[a];
	    	}
	    	I_ += NT;
	    	G_ += NU;
	    }
	}

	for (int j = 1; j <=  D; ++j) {
	    for (int i = 0; i <= C; ++i) {
		T * __restrict I_ = &I(0,i,j);
		U *G_ = &G(0,i);
#ifdef __IMTEL_COMPILER
#pragma ivdep
#pragma vector aligned
#endif
		for (int ij = 0; ij < M; ++ij) {
		    for(int a = 0; a < N; ++a) {
			G_[a] = dq*G_[a] + G_[a+MU];
			I_[a] = G_[a];
		    }
		    I_ += NT;
		    G_ += NU;
		}
		// for(int a = 0; a < MU; ++a) {
		//     G(a,i) = dq*G(a,i) + G(a,i+1);
		//     I(a,i,j) = G(a,i);
		// }
	    }

	    for (int i = C+1; i <= C+D-j; ++i) {
#ifdef __IMTEL_COMPILER
#pragma ivdep
#pragma vector aligned
#endif
		for(int a = 0; a < MU; ++a) {
		    G(a,i) = dq*G(a,i) + G(a,i+1);	    
		}
	    }
	}

#undef G
#undef I
    }

};


} // namespace kernel
} // namespace rysq

#endif /* RYSQ_KERNEL_TRANSFER_HPP_ */
