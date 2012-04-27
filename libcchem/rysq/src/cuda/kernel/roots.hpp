#ifndef _RYSQ_CUDA_KERNELS_ROOTS_H_
#define _RYSQ_CUDA_KERNELS_ROOTS_H_

#include "roots/rysq_roots.h"
#include "cuda/kernels/device.hpp"

namespace rysq {
    namespace cuda {
	namespace kernels {

	    template<ushort N> __device__
	    void roots(const Side<2,0> &bra, const Side<2,0> &ket,
		       const double2 *cAB, const double cutoff, const ushort Kkl,
		       double *t2W, double *C, const ushort &bRank) {

		using namespace device;

		/* parallel roots loop produces incorrect results */
		for (ushort Kij = 0; Kij < bra.K; ++Kij) {

		    double2 AB = cAB[Kij + Kkl*bra.K];
		    if (bRank == 0) C[Kij] = AB.x*bra.e[Kij]*ket.e[Kkl];
		    //if(fabs(C[Kij]) < cutoff) continue;
		    double AB1 = AB.y;

		    double *t2 = t2W + 2*N*Kij;
		    double *W = t2 + N;

		    if (N == 4) {
			    double rho = (bra.A[Kij])*(ket.A[Kkl])*AB1;
			    double X = rho*distance2(&bra.rA[Kij*3], &ket.rA[Kkl*3]);
			    rysq::roots4::evaluate(X, t2, W, bRank);
		    }
// 		    else 
// 			if (bRank == 0) {
// 			    double rho = (bra.A[Kij])*(ket.A[Kkl])*AB1;
// 			    double X = rho*distance2(&bra.rA[Kij*3], &ket.rA[Kkl*3]);
// 			    rysq::roots<N>(X, t2, W);
// 			}
 		    if(bRank < N) t2[bRank] *= AB1;

		}
		return;
	    }

	} // namespace kernels
    }
}

#endif /* _RYSQ_CUDA_KERNELS_ROOTS_H_ */
