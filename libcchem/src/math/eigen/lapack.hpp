#ifndef EIGEN_LAPACK_HPP
#define EIGEN_LAPACK_HPP

// #include <lapack.h>

namespace eigen {

    struct lapack {

	template<class M0, class V0, class M1>
	static void apply(M0 &A, V0 &w, M1 &V) {
	    char *jobz;
	    char *uplo;
	    integer *n;
	    doublereal *ap;
	    doublereal *w;
	    doublereal *z;
	    integer *ldz;
	    doublereal *work;
	    integer *info;

	    dspev_(&jobz, &uplo, &n, ap, w, z, &ldz, work, &info);

    };


}

#endif // EIGEN_LAPACK_HPP
