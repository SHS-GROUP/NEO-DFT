#include "rysq-types.hpp"

#include "cuda/kernel/kernel.hpp"
#include "cuda/kernel/eri.hpp"
#include "cuda/kernel/fock.hpp"

#include "cuda/kernel/quadrature.hpp"
#include "cuda/kernel/quadrature2.hpp"

#include <boost/preprocessor/seq/for_each_product.hpp>
#include <boost/preprocessor/seq/elem.hpp>
#include <boost/preprocessor/seq/enum.hpp>

namespace rysq {
namespace cuda {
namespace kernel {


template<rysq::type A, rysq::type B, rysq::type C, rysq::type D, class Transform>
kernel::Eri<Transform>* instantiate(const Shell::Quartet &quartet,
				    const Context &context) {
    int N = quartet.L()/2 + 1;
    try {
    	{//if (quartet.size() <= 10) {// || N*quartet.K() >= 512) {
    	    return new Kernel2<A, B, C, D, Transform>(quartet, context);
    	}
    }
    catch (kernel::kernel_not_implemented) {}
    try {
    	return new Kernel1<A, B, C, D, Transform>(quartet, context);
    }
    catch (kernel::kernel_not_implemented) {}
    return NULL;
}


#define ELEM BOOST_PP_SEQ_ELEM
#define ENUM BOOST_PP_SEQ_ENUM
#define FOR_EACH_PRODUCT BOOST_PP_SEQ_FOR_EACH_PRODUCT
#define KERNEL(r, params)					\
    template<>							\
    kernel::Eri<ELEM(4, params)>*				\
    new_<ENUM(params)>(const Shell::Quartet &quartet,		\
		       const Context &context) {		\
	return instantiate<ENUM(params)>(quartet, context);	\
    }

    // generate AB kernels
    FOR_EACH_PRODUCT(KERNEL, \
		     (RYSQ_CUDA_KERNEL_SHELL_A) \
		     (RYSQ_CUDA_KERNEL_SHELL_B) \
		     (RYSQ_TYPES)(RYSQ_TYPES)   \
		     ((eri::Transform)(fock::Transform<>)));


#undef KERNEL
#undef FOR_EACH_PRODUCT
#undef ELEM
#undef ENUM

} // namespace kernel
} // namespace cuda
} // namespace rysq

