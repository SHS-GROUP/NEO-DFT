#ifndef RYSQ_KERNEL_ERI2_HPP
#define RYSQ_KERNEL_ERI2_HPP

#include "kernel/eri.hpp"
#include "kernel/quadrature2.hpp"
#include "kernel/vector.hpp"
#include "transpose.hpp"

#include <boost/utility/enable_if.hpp>

namespace rysq {
namespace kernel {

template<type T0, type T1, type T2, type T3>
struct Eri<meta::state<T0, T1>, meta::state<T2, T3>,
	   typename boost::enable_if<
	       //boost::mpl::bool_<meta::sorted<T0,T1,T2,T3>::braket::L == 0>
	       quadrature::impl<
	        	   typename meta::sorted<T0,T1,T2,T3>::braket>
	       >::type
	   >
    : public Eri<>
{
    static const bool value = true;
    typedef meta::state<T0, T1> bra;
    typedef meta::state<T2, T3> ket;
    typedef typename meta::sorted<T0,T1,T2,T3>::braket braket;
    static const int transpose_value = meta::sorted<T0,T1,T2,T3>::value;

    typedef kernel::Transform<bra,ket> Transform;
    typedef Eri<>::Data Data;

    Eri(const Quartet<Shell> &quartet, Transform *transform)
	: Eri<>(transpose(quartet, transpose_value)),
	  transform_(transform) {
	// std::cout << quartet << transpose(quartet, transpose_value) << std::endl;	
    }

    ~Eri() { delete transform_; }

    void operator()(const Quartet<Center> &r,
		    Data &data,
		    const rysq::Eri::Parameters &parameters) {
	const Quartet<Center> r_ = transpose(r, transpose_value);
	double scale = rysq::SQRT_4PI5;
	double cutoff = parameters.cutoff/(scale*this->quartet_.K());
	cutoff *= CUTOFF_SCALE;
	double Q[braket::size]  __attribute__ ((aligned(16))) = { 0.0 };
	// std::cout << quartet_ << transpose_value<< std::endl;	
	typedef kernel::vector<3> vector;
	quadrature::apply<braket>(this->quartet_,
				  vector(r_[0]), vector(r_[1]),
				  vector(r_[2]), vector(r_[3]),
				  Q, transpose_value, cutoff);
	((*transform_)(data))(Q, scale);
    }

private:
    Transform *transform_;
};


} // namespace kernel
} // namespace rysq

#endif /* RYSQ_KERNEL_ERI2_HPP */
