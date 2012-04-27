#ifndef RYSQ_CUDA_KERNEL_KERNEL_HPP
#define RYSQ_CUDA_KERNEL_KERNEL_HPP

#include <rysq/cuda.hpp>
#include "cuda/detail.hpp"
#include "boost/cuda/device/device.hpp"

namespace rysq {
namespace cuda {
namespace kernel {

typedef boost::cuda::device::Block thread_block;

template<class Transform>
struct Eri {
    static Eri* new_(const detail::Quartet &quartet,
		     const rysq::Transpose &transpose);
    virtual ~Eri() {}
    virtual void operator()(const detail::Centers &centers,
			    const detail::Quartets &quartets,
			    const cuda::Eri::Parameters &p,
			    const boost::cuda::stream &stream,
			    Transform transform,
			    ushort threads = 0, size_t shared = 0) = 0;
};

struct Transform;

template<class, class Transform, typename T = double, class Enable = void>
struct quadrature {
    static const bool value = false;
};


template<class Parameters, class Transform, class Enable = void>
struct Kernel {
    static Eri<Transform>* new_(const detail::Quartet &quartet,
				const rysq::Transpose &transpose) {
	throw std::runtime_error("kernel specialization not found");
    }
};

} // namespace kernel
} // namespace cuda
} // namespace rysq


#include "cuda/kernel/quadrature.hpp"
#include "cuda/kernel/quadrature2.hpp"

#include <boost/preprocessor/seq/for_each_product.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/seq/enum.hpp>
#include <boost/preprocessor/seq/elem.hpp>

namespace rysq {
namespace cuda {
namespace kernel {

template<class Transform>
kernel::Eri<Transform>*
kernel::Eri<Transform>::new_(const detail::Quartet &quartet,
			     const rysq::Transpose &transpose) {

#define ELEM BOOST_PP_SEQ_ELEM
#define ENUM BOOST_PP_SEQ_ENUM
#define FOR_EACH_PRODUCT BOOST_PP_SEQ_FOR_EACH_PRODUCT

#define KERNEL(r, params)						\
    if ((quartet[0] == ELEM(0, params)) &&				\
	(quartet[1] == ELEM(1, params)) &&				\
	(quartet[2] == ELEM(2, params)) &&				\
	(quartet[3] == ELEM(3, params))) {				\
	try {								\
	    typedef meta::braket<ENUM(params)> braket;			\
	    return Kernel<braket, Transform>::new_(quartet, transpose);	\
	}								\
	catch (std::exception&) { }					\
    } else
    //FOR_EACH_PRODUCT(KERNEL, (RYSQ_TYPES)(RYSQ_TYPES)(RYSQ_TYPES)(RYSQ_TYPES)) {}
#undef KERNEL

    if (quartet.size() > 2400) {
	//std::cout << "NULL" << std::endl;
	return NULL;
    }
    if (quartet.size() < 16) {
	//std::cout << "NULL" << std::endl;
	return NULL;
    }

    int N = quartet.L()/2 + 1;
    int mask = cxx::utility::bitmask(quartet[0].L, quartet[1].L,
				     quartet[2].L, quartet[3].L);

#define SEQ_MASK        (Q1111)(Q1110)(Q1011)(Q1100)(Q1010)
#define SEQ_N		(2)(3)(4)(5)

#define KERNEL(r, params)						\
    if ((mask == ELEM(0, params)) && (N == ELEM(1, params))) {		\
	try {								\
	    return Kernel<Parameters<ENUM(params)>,			\
		Transform>::new_(quartet, transpose);			\
	}								\
	catch (std::exception&) { return NULL; }			\
    } else

    FOR_EACH_PRODUCT(KERNEL, (SEQ_MASK)(SEQ_N)) { }

#undef SEQ_N
#undef KERNEL
#undef SEQ_MASK

#undef ELEM
#undef ENUM
#undef FOR_EACH_PRODUCT

    {
	//std::cout << "quadrature " << N << " " << mask << std::endl;
	return NULL;
    }

}

} // namespace kernel
} // namespace cuda
} // namespace rysq


#endif // RYSQ_CUDA_KERNEL_KERNEL_HPP
