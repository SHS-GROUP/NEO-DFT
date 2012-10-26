#include "cuda/kernel/kernel.hpp"
#include "cuda/kernel/fock.hpp"

#include "boost/utility/profiler.hpp"

namespace rysq {
namespace cuda {
namespace detail {

    typedef kernel::fock::Transform<> Transform;
    typedef kernel::Eri<Transform> Kernel;

    void Fock::operator()(const Centers &centers,
			  const Quartets &quartets,
			  const Set &set, Mutex mutex,
			  const cuda::Fock::Parameters &p,
			  cudaStream_t stream) {

	// BOOST_PROFILE_LINE;

	double scale[] = { p.cscale, p.xscale };
	Transform T(set, scale, mutex);

	Kernel *impl = static_cast<Kernel*>(impl_);
	(*impl)(centers, quartets, cuda::Eri::Parameters(p.cutoff), stream, T);
    }


    Fock::Fock(const detail::Shell::Quartet &quartet,
	       const Context &context,
	       const Transpose &transpose) 
	: quartet_(quartet)
    {
	BOOST_PROFILE_LINE;
	//std::cout << quartet << std::endl;
	this->impl_ = kernel::new_<Transform>(this->quartet_, context);
	if (!this->impl_)
	    throw rysq::cuda::exception("no kernel");
    }

    Fock::~Fock() {
	BOOST_PROFILE_LINE;
	delete static_cast<Kernel*>(this->impl_);
    }

}
}
}
