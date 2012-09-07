#include "kernel/eri.hpp"

#include "rysq-core.hpp"

#include "kernel/eri1.hpp"
#include "kernel/eri2.hpp"
#include "kernel/quadrature1.hpp"
#include "kernel/quadrature2.hpp"

#include "kernel/quadrature2-impl.hpp"
#include "kernel/new.hpp"


namespace rysq {
namespace kernel {


#define TYPES	(rysq::SP)(rysq::S)(rysq::P)(rysq::D)(rysq::F)

#define ERI(r, types) 							\
    kernel::find<BOOST_PP_SEQ_ENUM(types)>::type; 

    // BOOST_PP_SEQ_FOR_EACH_PRODUCT(ERI, (TYPES)(TYPES)(TYPES)(TYPES))

#undef ERI
#undef TYPES

namespace instantiate_eri {

    template<class bra, class ket>
    struct Transform : kernel::Transform<bra,ket> {
	typedef kernel::Transform<bra,ket> Base;
	Base& operator()(typename Base::Data &data) { return *this; }
	void operator()(const double *Q, double scale) {}
    };

    template<class bra>
    struct Transform <bra, void> : kernel::Transform<bra,void> {
	typedef kernel::Transform<bra> Base;
	Base& operator()(typename Base::Data &data) { return *this; }
	void operator()(int k,int l,int kl,
			const double *Q, double scale) {}
    };

    void instance( const Quartet < Shell> &quartet) {
	delete  kernel::new_< Transform>(quartet);

    }
}


} // namespace kernel
} // namespace rysq

