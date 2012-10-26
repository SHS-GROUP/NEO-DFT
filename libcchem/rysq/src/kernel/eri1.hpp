#ifndef RYSQ_KERNEL_ERI1_HPP
#define RYSQ_KERNEL_ERI1_HPP

#include "kernel/quadrature1.hpp"
#include <math.h>
#include <boost/mpl/int.hpp>
#include <boost/mpl/assert.hpp>

#include "rysq-core.hpp"
#include "kernel/eri.hpp"
#include "vector.hpp"


namespace rysq {
namespace kernel {


template<typename T, size_t N>
struct align {
#define ALIGN_(n) ((A/sizeof(T) - (n)%(A/sizeof(T)))%(A/sizeof(T)))
    static const size_t A = 16;
    static const size_t value = ALIGN_(N);
    static size_t get(size_t M) { return ALIGN_(M); }
#undef ALIGN_
};

template<class bra_, int N>
struct Eri<bra_, boost::mpl::int_<N> > : public Eri<>
{
    //BOOST_MPL_ASSERT_MSG((bra_::L > 0), );
    typedef bra_ bra;
    typedef void ket;
    typedef kernel::Transform<bra> Transform;
    typedef Eri<>::Data Data;
    Eri(const Quartet<Shell> &quartet, Transform *transform)
	: Eri<>(quartet), transform_(transform)
    {
	primitives_.allocate<align>(quartet);
    }

    ~Eri() { delete transform_; }

    virtual void operator()(const Quartet<Center> &r, Data &data,
			    const rysq::Eri::Parameters &parameters) {
	double scale = rysq::SQRT_4PI5;
	double cutoff = parameters.cutoff/(scale*this->quartet_.K());
	cutoff *= CUTOFF_SCALE;
	typedef kernel::vector<3> vector;
	quadrature::apply<bra,N,align>(this->quartet_,
				       vector(r[0]), vector(r[1]),
				       vector(r[2]), vector(r[3]),
				       scale, cutoff, primitives_,
				       (*transform_)(data));
    }

private:
    Transform *transform_;
    //quadrature::Primitives<long double, long double> primitives_; 
    quadrature::Primitives<double, double> primitives_;
    //quadrature::Primitives<float, double> primitives_;
    //  quadrature::Primitives<float, float> primitives_;

};


} // namespace kernel
} // namespace rysq

#endif // RYSQ_KERNEL_ERI1_HPP

