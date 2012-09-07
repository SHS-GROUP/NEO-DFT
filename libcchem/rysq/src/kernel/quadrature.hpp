#ifndef RYSQ_KERNEL_QUADRATURE_HPP_
#define RYSQ_KERNEL_QUADRATURE_HPP_

#include "rysq-core.hpp"
#include "vector.hpp"


#include <memory>
#include <assert.h>
#include <stdlib.h>
#include <math.h>


namespace rysq {
namespace kernel {
namespace quadrature {


template<int N>
struct roots {
    static const size_t value = N;
};

template<class braket>
void apply(const Quartet<Shell> &quartet,
	   const vector<3> &ri, const vector<3> &rj,
	   const vector<3> &rk, const vector<3> &rl,
	   double (&Q)[braket::size], int transpose,
	   double cutoff);

template<typename T, typename U>
struct Primitives {
    int Kmax, K;
    int size2d;
    void *memory16_;
    U *Int2d1, *Gx_, *Gy_, *Gz_;
    T *Ints2d;
    double *C[4];

    template<template<typename, size_t> class align>
    void allocate(const Quartet<Shell> &quartet);

    ~Primitives() { if (memory16_) free(memory16_); }

    U* Gx(int K) { return (Gx_) ? Gx_ : get<U,0>(K); }
    U* Gy(int K) { return (Gy_) ? Gy_ : get<U,1>(K); }
    U* Gz(int K) { return (Gz_) ? Gz_ : get<U,2>(K); }

    U* transfer() { return Int2d1; }

     T* Ix(int K = 0) { return get<T,0>(K); }
     T* Iy(int K = 0) { return get<T,1>(K); }
     T* Iz(int K = 0) { return get<T,2>(K); }

     const T* Ix(int K = 0) const { return get<T,0>(K); }
     const T* Iy(int K = 0) const { return get<T,1>(K); }
     const T* Iz(int K = 0) const { return get<T,2>(K); }

private:
    template<typename R, size_t q> R* get(int i) {
	return reinterpret_cast<R*>(&Ints2d[(3*i + q)*size2d]);
    }
    template<typename R, size_t q>  const R* get(int i)const {
	return (const_cast<Primitives*>(this))->get<R,q>(i);
    }
};


template<class bra, size_t N, template<typename, size_t> class align,
	 class Transform, typename T, typename U>
 void apply(const Quartet<Shell> &quartet,
	    const vector<3> &ri, const vector<3> &rj,
	    const vector<3> &rk, const vector<3> &rl,
	    double scale, double cutoff,
	    Primitives<T,U> &primitives,
	    Transform &transform);

template<typename T, typename U>
template< template<typename, size_t> class align>
void Primitives<T,U>::allocate(const Quartet<Shell> &quartet) {
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

    size_t NU = N + align<U,0>::get(N);
    int size2d0 = NU*(a.L + b.L + 1)*(c.L + d.L + 1);
    int size2d1 = NU*(a.L + 1)*(b.L + 1)*(c.L + d.L + 1);

    size_t NT = N + align<T,0>::get(N);
    size2d = NT*(a.L + 1)*(b.L + 1)*(c.L + 1)*(d.L + 1);

    //int nc = a->nc*b->nc*c->nc*d->nc;
    //int K = a->K*b->K*c->K*d->K;

    int nc = quartet.nc();
    int K = quartet.K();

    Kmax = K;//(RYSQ_WITH_INT2D_BLOCK*1024)/((3*size2d+nc)*sizeof(double));
    Kmax = std::max(Kmax, 1);
    Kmax = std::min(Kmax, K);

    bool reuse = ((transfer != 1) && (size2d*sizeof(T) >= size2d0*sizeof(U)));
    size_t size_tmp = 0;
    size_tmp += (!reuse) ? 3*size2d0 : 0;
    size_tmp += (transfer == 2) ? size2d1 : 0;

    size_t size = 0;
    size += (3*size2d)*Kmax;
    size_t bytes = size_tmp*sizeof(U) + size*sizeof(T) + nc*Kmax*sizeof(double);
    int status = posix_memalign(&memory16_, 16, bytes);

    assert(status == 0);
    if (status != 0) throw std::bad_alloc();

    U *ptr = (U*)memory16_;
    if (!reuse) {
	Gx_ = ptr + 0*size2d0;
	Gy_ = ptr + 1*size2d0;
	Gz_ = ptr + 2*size2d0;
	ptr += 3*size2d0;
    }

    if ( transfer == 2) {
    Int2d1 = ptr;
    ptr += size2d1;
    }

    T *ptr64 = (T*)ptr;
    Ints2d = ptr64;
    ptr64 += 3*size2d*Kmax;

    for (int i = 0; i < (c.nc*d.nc); ++i) {
	C[i] = (double*)ptr64 + i*(a.nc*b.nc)*Kmax;
    }
}


} // namespace rysq {
} // namespace quadrature
} // namespace kernel


#endif /* RYSQ_KERNEL_QUADRATURE_HPP_ */
