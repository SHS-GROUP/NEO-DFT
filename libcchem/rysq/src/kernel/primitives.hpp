#ifndef _RYSQ_KERNEL_PRIMITIVES_HPP_
#define _RYSQ_KERNEL_PRIMITIVES_HPP_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <memory>
#include <assert.h>
#include <stdlib.h>
#include <math.h>

#include <rysq/core.hpp>

#define pad2(N) ((N)+(N)%2)

BEGIN_NAMESPACE(rysq, kernel, eri_)

template<typename R>
struct Primitives {
    int Kmax, K;
    bool final;
    int size2d;
    void *memory16_;
    R *Int2d1, *Gx_, *Gy_, *Gz_;
    double *Ints2d;
    double *C[4];

    Primitives(const Quartet<Shell> &quartet);
    ~Primitives() { if (memory16_) free(memory16_); }

    template<typename T> T* Gx(int K) { return (Gx_) ? (T*)Gx_ : (T*)Ix(K); }
    template<typename T> T* Gy(int K) { return (Gy_) ? (T*)Gy_ : (T*)Iy(K); }
    template<typename T> T* Gz(int K) { return (Gz_) ? (T*)Gz_ : (T*)Iz(K); }

    template<typename T> T* transfer() { return (T*)Int2d1; }

    double* Ix(int K = 0) { return &Ints2d[(3*K + 0)*size2d]; }
    double* Iy(int K = 0) { return &Ints2d[(3*K + 1)*size2d]; }
    double* Iz(int K = 0) { return &Ints2d[(3*K + 2)*size2d]; }

    const double* Ix(int K = 0) const { return &Ints2d[(3*K + 0)*size2d]; }
    const double* Iy(int K = 0) const { return &Ints2d[(3*K + 1)*size2d]; }
    const double* Iz(int K = 0) const { return &Ints2d[(3*K + 2)*size2d]; }

};

template<typename R>
inline Primitives<R>::Primitives(const Quartet<Shell> &quartet) {
    memory16_ = NULL;
    Gx_ = Gy_ = Gz_ = NULL;
    Int2d1 = NULL;
    Ints2d = NULL;
    C[0] = C[1] = C[2] = C[3] = NULL;

    // assert(quartet.L() > 0);
    int N = quartet.L()/2 + 1;

    const Shell &a = quartet[0], &b = quartet[1], &c = quartet[2], &d = quartet[3];
    // quartet.unpack(a, b, c, d);

    int transfer = (a.L && b.L) + (c.L && d.L);

    int size2d0 = pad2(N)*(a.L + b.L + 1)*(c.L + d.L + 1);
    int size2d1 = pad2(N)*(a.L + 1)*(b.L + 1)*(c.L + d.L + 1);

    size2d = pad2(N)*(a.L + 1)*(b.L + 1)*(c.L + 1)*(d.L + 1);

    //int nc = a->nc*b->nc*c->nc*d->nc;
    //int K = a->K*b->K*c->K*d->K;

    int nc = quartet.nc();
    int K = quartet.K();

    Kmax = (RYSQ_WITH_INT2D_BLOCK*1024)/((3*size2d+nc)*sizeof(double));
    Kmax = std::max(Kmax, 1);
    Kmax = std::min(Kmax, K);

    size_t size_tmp = 0;
    size_tmp += (transfer == 1) ? 3*size2d0 : 0;
    size_tmp += (transfer == 2) ? size2d1 : 0;

    size_t size = 0;
    size += (3*size2d + nc)*Kmax;
    size_t bytes = size_tmp*sizeof(R) + size*sizeof(double);
    // std::cout << K << " "  << nc << "\n";
    assert(posix_memalign(&memory16_, 16, bytes) == 0);

    R *ptr = (R*)memory16_;
    if (transfer == 1) {
	Gx_ = ptr + 0*size2d0;
	Gy_ = ptr + 1*size2d0;
	Gz_ = ptr + 2*size2d0;
	ptr += 3*size2d0;
    }
    else if (transfer == 2) {
	Int2d1 = ptr;
	ptr += size2d1;
    }

    double *ptr64 = (double*)ptr;
    Ints2d = ptr64;
    ptr64 += 3*size2d*Kmax;

    for (int i = 0; i < (c.nc*d.nc); ++i) {
	C[i] = ptr64 + i*(a.nc*b.nc)*Kmax;
    }
}

END_NAMESPACE(rysq, kernel, eri_)

#endif // _RYSQ_KERNEL_PRIMITIVES_HPP_

