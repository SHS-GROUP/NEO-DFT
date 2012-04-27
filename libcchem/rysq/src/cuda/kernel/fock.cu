#define KERNEL_GLOBAL_PREFIX fock
// #define RYSQ_CUDA_KERNEL_USE_CONSTANT 1

#include "cuda/kernel/kernel.hpp"
#include "cuda/kernel/fock.hpp"
#include <boost/preprocessor/repetition/repeat.hpp>

namespace rysq {
namespace cuda {
namespace detail {

    typedef kernel::fock::Transform<> Transform;
#define KERNEL kernel::Eri<Transform>

    void Fock::operator()(const Centers &centers,
			  const Quartets &quartets,
			  const Set &set, int* const (&mutex)[6],
			  const cuda::Fock::Parameters &p,
			  const boost::cuda::stream &stream) {

	double scale[] = { p.cscale, p.xscale };
	Transform T(set, scale, mutex);

	ushort threads = 0;
	using rysq::detail::index;
	
#define MAX(z, i, data)							\
	threads = std::max<size_t>(threads, T.matrix.size(fock_index::ij<i>::type()));
	BOOST_PP_REPEAT(6, MAX, ())
#undef MAX

	size_t size = 1;
	for (int i = 0; i < 4; ++i) {
	    size *= T.matrix.block[i];
	}

	// std::cout <<  shared<<std::endl;
	size_t shared = (size + 2*threads)*sizeof(double);

	KERNEL *impl = static_cast<KERNEL*>(impl_);
	(*impl)(centers, quartets, cuda::Eri::Parameters(p.cutoff),
		stream, T, threads, shared);
    }


    Fock::Fock(const rysq::Quartet<rysq::Shell> &quartet,
	       const Transpose &transpose) 
	: quartet_(quartet)
    {
	//std::cout << quartet << std::endl;
	this->impl_ = KERNEL::new_(this->quartet_, transpose);
	if (!this->impl_) throw std::exception();
    
    }

    Fock::~Fock() {
	delete static_cast<KERNEL*>(this->impl_);
    }

}
}
}
