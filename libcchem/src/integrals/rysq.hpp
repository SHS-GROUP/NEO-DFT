#ifndef INTEGRALS_RYSQ_HPP
#define INTEGRALS_RYSQ_HPP

#include "basis/basis.hpp"
#include "integrals/quartet.hpp"
#include "exception.hpp"

#include <rysq.hpp>
#include "adapter/rysq.hpp"

#ifdef RYSQ_CUDA
#include "cuda.hpp"
#define CCHEM_INTEGRALS_ERI_CUDA
#endif

#include <boost/array.hpp>
#include <boost/ptr_container/ptr_map.hpp>
#include "boost/utility/profiler.hpp"

namespace integrals {
namespace rysq {

    inline void initialize() {
	::rysq::initialize();
    }

    struct Eri {
	struct Cuda;
	Eri(Quartet<const Basis::Shell::Data&> quartet) :
	    A_(quartet.p), B_(quartet.q),
	    C_(quartet.r), D_(quartet.s),
	    eri_(::rysq::Quartet< ::rysq::Shell >(A_, B_, C_, D_)) {}
	void operator()(const std::vector<Basis::Center> &centers,
			const std::vector< boost::array<int,4> > &quartets,
			double *G, double cutoff) {
	    //BOOST_PROFILE_FUNCTION();
	    // std::cout <<  "integralss: " << quartets.size() << std::endl;
	    eri_(centers, quartets, G, cutoff);
	}
    private:
	adapter::rysq::Shell A_, B_, C_, D_;
	::rysq::Eri eri_;
    };

#ifdef RYSQ_CUDA
    struct Eri::Cuda {
	struct not_implemented : std::runtime_error {
	    not_implemented(const char *what) : std::runtime_error(what) {}
	};
	Cuda(const Basis &basis)
	    : context_()
	{
	    std::vector< ::rysq::Center > centers;
	    foreach (const Basis::Center &c, basis.centers()) {
		centers.push_back(c);
	    }
	    this->centers_.assign(centers);
	}
	void operator()(const Basis::Shell::Data &A,
			const Basis::Shell::Data &B,
			const Basis::Shell::Data &C,
			const Basis::Shell::Data &D,
			const std::vector< boost::array<int,4> > &quartets,
			double *I, double cutoff,
			cuda::stream stream) {
	    BOOST_PROFILE_LINE;
	    ::rysq::cuda::Eri eri_
		(this->insert(A),
		 this->insert(B),
		 this->insert(C),
		 this->insert(D),
		 this->context_);
	    if (!eri_) throw CCHEM_EXCEPTION("no ::rysq::cuda::Eri kernel");
	    ::rysq::cuda::Stream s(stream);
	    this->quartets_.assign(quartets, s);
	    ::rysq::cuda::Eri::Parameters p;
	    p.cutoff = cutoff;
	    eri_(this->centers_, this->quartets_, I, p, s);
	}
    private:
	::rysq::cuda::Context context_;	    
	::rysq::cuda::Centers centers_;
	boost::ptr_map<Basis::Shell::Data::Key, ::rysq::cuda::Shell> shells_;
	::rysq::cuda::Quartets quartets_;
	::rysq::cuda::Shell& insert(const Basis::Shell::Data &S) {
	    Basis::Shell::Data::Key k = S.key();
	    if (!shells_.count(k)) {
		shells_.insert(k, new ::rysq::cuda::Shell(adapter::rysq::Shell(S)));
	    }
	    assert(shells_.find(k) != shells_.end());
	    assert(shells_.find(k)->second);
	    return *(shells_.find(k)->second);
	}
    };
#endif // RYSQ_CUDA

} // rysq
} // integrals


#endif // INTEGRALS_RYSQ_HPP
