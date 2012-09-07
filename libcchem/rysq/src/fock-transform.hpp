#ifndef RYSQ_FOCK_TRANSFORM_HPP_
#define RYSQ_FOCK_TRANSFORM_HPP_


namespace rysq {
namespace fock {


template<class bra = void, class ket = void>
struct Transform;

template<>
struct Transform<void, void> {
    struct Data : kernel::Transform<>::Data {
	typedef boost::array<const double*,6> density_type;
	typedef boost::array<double*,6> fock_type;
	Data() {}
	Data(density_type density, fock_type fock)
	    : density(density), fock(fock) {}
	density_type density; fock_type fock;
    };
};

typedef Transform<>::Data Data;

}} // namespace

#include "kernel/eri.hpp"
#include "fock-eval.hpp"


namespace rysq {
namespace fock {


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
	eval<ni, nj, nk, nl>(data_.density, data_.fock, Q, scale);
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
	// std::cout << scale << std::endl;
	const size_t ni = bra::A::size;
	const size_t nj = bra::B::size;
	// std::cout << ni << " " << nj << std::endl;
	eval<ni, nj>(k, l, kl, data_.density, data_.fock, Q, scale);
    }
private:
    Data data_;
};


}} // namespace


#endif /* RYSQ_FOCK_TRANSFORM_HPP_ */
