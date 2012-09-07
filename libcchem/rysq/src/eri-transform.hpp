#ifndef RYSQ_ERI_TRANSFORM_HPP_
#define RYSQ_ERI_TRANSFORM_HPP_

#include "kernel/eri.hpp"

namespace rysq {
namespace eri {

template<class bra = void, class ket = void>
struct Transform;

template<>
struct Transform<void, void> {
    struct Data : kernel::Transform<>::Data {
	Data() {}
	Data(double *Q) : data_(Q) {}
	double* begin() { return data_; }
    private:
	double *data_;
    };
};

template<class bra, class ket>
struct Transform : kernel::Transform<bra,ket> {
    typedef kernel::Transform<bra,ket> Base;
    typedef Transform<>::Data Data;
    Base& operator()(typename Base::Data &data) {
	this->data_ = static_cast<Data&>(data);
	return *this;
    }
    void operator()(const double *Q, double scale) {
	const size_t ni = bra::A::size, nj = bra::B::size;
	const size_t nk = ket::A::size, nl = ket::B::size;
	const size_t size = ni*nj*nk*nl;
	double *data = data_.begin();
	for (size_t i = 0; i < size; ++i) {
	    data[i] += scale*Q[i];
	}
    }
private:
    Data data_;
};


template<class bra>
struct Transform<bra,void> : kernel::Transform<bra> {
    typedef kernel::Transform<bra> Base;
    typedef Transform<>::Data Data;
    Base& operator()(typename Base::Data &data) {
	this->data_ = static_cast<Data&>(data);
	return *this;
    }
    void operator()(int k,int l,int kl, const double *Q, double scale) {
	const size_t ni = bra::A::size;
	const size_t nj = bra::B::size;
	const size_t size = ni*nj;
	double *data = data_.begin() + kl*size; 
	for (size_t i = 0; i < size; ++i) {
	    data[i] += scale*Q[i];
	}
    }
private:
    Data data_;
};

}} // namespace

#endif /* RYSQ_ERI_TRANSFORM_HPP_ */
