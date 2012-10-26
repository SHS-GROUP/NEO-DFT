#include <rysq.hpp>
#include "meta.hpp"

using namespace rysq;

template<class  M>
class Transform {
    
    M C1_, C2_;
    M T2_;

    template<rysq::type A, rysq::type B,size_t N>
    void operator()(int k, int l, const double (&Q)[N]) {
	const size_t ni = meta::shell<A>::size;
	const size_t nj = meta::shell<B>::size;
	size_t size1,size2;
	for (int b = 0; b < size1; ++b) {
	    double *C1 = &C1_.block(0,0);
	    double T1[nj];
	    for (int j = 0; j < nj; ++j) {
	        for (int i = 0; i < ni; ++i) {
		    T1[j] = C1[i]*Q[i+j*ni];
		}
	    }

	    double *C2 = &C2_.data()[0];
	    double *T2 = &T2_.iterator2()[a];
	    for (int a = 0; a < size2; ++a) {
	        for (int j = 0; j < nj; ++j) {
		    T2[a] = C2[j+a*nj]*T1[j];
		}
	    }
	}
    }

};
