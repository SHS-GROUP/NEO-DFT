#ifndef RYSQ_CUDA_KERNEL_KERNEL_HPP
#define RYSQ_CUDA_KERNEL_KERNEL_HPP

#include "cuda/detail.hpp"
#include "cuda/kernel/quartet.hpp"

#include <boost/preprocessor/seq/for_each_product.hpp>
#include <boost/preprocessor/seq/elem.hpp>
#include <boost/preprocessor/seq/enum.hpp>


namespace rysq {
namespace cuda {
namespace kernel {

    struct kernel_not_implemented {};

    template<class Transform>
    struct Eri {
	virtual ~Eri() {}
	virtual void operator()(const detail::Centers &centers,
				const detail::Quartets &quartets,
				const cuda::Eri::Parameters &p,
				cudaStream_t stream,
				Transform transform) {
	    throw kernel_not_implemented();
	}
    };

    struct Dim3 : dim3 {
	explicit Dim3(int x = 1, int y = 1, int z = 1) : dim3(x,y,z) {}
	int size() const { return dim3::x*dim3::y*dim3::z; }
    };

    inline int align(int value, int size) {
	while (value%size) ++value;
	return value;
    }


template<rysq::type A, rysq::type B, rysq::type C, rysq::type D, class Transform>
kernel::Eri<Transform>* new_(const Shell::Quartet &quartet,
			     const Context &context);

template<class Transform>
kernel::Eri<Transform>* new_(Shell::Quartet quartet,
			     const Context &context) {
 
    quartet.sort();

#define ELEM BOOST_PP_SEQ_ELEM
#define ENUM BOOST_PP_SEQ_ENUM
#define FOR_EACH_PRODUCT BOOST_PP_SEQ_FOR_EACH_PRODUCT
#define TYPES (rysq::SP)(rysq::S)(rysq::P)(rysq::D)(rysq::F)

#define KERNEL(r, params)					\
    if ((quartet[0] == ELEM(0, params)) &&			\
	(quartet[1] == ELEM(1, params)) &&			\
	(quartet[2] == ELEM(2, params)) &&			\
	(quartet[3] == ELEM(3, params))) {			\
	return new_<ENUM(params), Transform>(quartet, context);	\
    }								\
    else

    // generate kernels
    FOR_EACH_PRODUCT(KERNEL, (TYPES)(TYPES)(TYPES)(TYPES)) {}

#undef KERNEL
#undef TYPES
#undef FOR_EACH_PRODUCT
#undef ENUM
#undef ELEM

    {
	return NULL;
    }

}


} // namespace kernel
} // namespace cuda
} // namespace rysq


#endif // RYSQ_CUDA_KERNEL_KERNEL_HPP
