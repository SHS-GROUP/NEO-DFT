#include "cuda/kernel/kernel.hpp"
#include "cuda/kernel/eri.hpp"

namespace rysq {
namespace cuda {

    typedef kernel::eri::Transform Transform;
    typedef kernel::Eri<Transform> Kernel;

    void detail::Eri::operator()(const detail::Centers &centers,
    				 const detail::Quartets &quartets,
    				 List I,
    				 const cuda::Eri::Parameters &p,
    				 cudaStream_t stream) {
    	typedef kernel::eri::Transform T;
    	Kernel *impl = static_cast<Kernel*>(impl_);
    	(*impl)(centers, quartets, p, stream, T(I));
    }

    detail::Eri::Eri(const detail::Shell::Quartet &quartet,
		     const Context &context,
    		     const rysq::Transpose &transpose)
    	: quartet_(quartet)
    {
    	//std::cout << quartet << std::endl;
    	this->impl_ = kernel::new_<Transform>(this->quartet_, context);
    	if (!this->impl_) throw std::exception();
    }

    detail::Eri::~Eri() {
    	delete static_cast<Kernel*>(this->impl_);
    }

}
}
