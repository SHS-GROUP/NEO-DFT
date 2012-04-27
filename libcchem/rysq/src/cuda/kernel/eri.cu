#define KERNEL_GLOBAL_PREFIX eri

#include "cuda/kernel/kernel.hpp"
#include "cuda/kernel/eri.hpp"

namespace rysq {
namespace cuda {
namespace detail {

#define KERNEL kernel::Eri<kernel::eri::Transform>

    void Eri::operator()(const detail::Centers &centers,
			 const detail::Quartets &quartets,
			 List I,
			 const cuda::Eri::Parameters &p,
			 const boost::cuda::stream &stream) {
	typedef kernel::eri::Transform T;
	KERNEL *impl = static_cast<KERNEL*>(impl_);
	(*impl)(centers, quartets, p, stream, T(I));
    }

    Eri::Eri(const rysq::Quartet<rysq::Shell> &quartet,
	     const rysq::Transpose &transpose)
	: quartet_(quartet)
    {
	//std::cout << quartet << std::endl;
	this->impl_ = KERNEL::new_(this->quartet_, transpose);
	if (!this->impl_) throw std::exception();
    }

    Eri::~Eri() {
	delete static_cast<KERNEL*>(this->impl_);
    }

}
}
}
