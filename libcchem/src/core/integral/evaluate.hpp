#ifndef INTEGRAL_EVALUATE_HPP
#define INTEGRAL_EVALUATE_HPP

#include "basis/basis.hpp"
#include "core/integral/aux.hpp"

#include <rysq.hpp>
#include "adapter/rysq.hpp"
#ifdef RYSQ_CUDA
#define HAVE_GPU_INTEGRALS
#endif

#include <vector>
#include <boost/array.hpp>
#include "boost/utility/profiler.hpp"

namespace integral {
  
    struct rysq_ {

	typedef Basis::Shell Shell;
	//typedef Basis::Shell::Data UniqueShell;
	typedef Basis::Center Center;

	static void evaluate(const Shell &A, const Shell &B,
			     const Shell &C, const Shell &D,
			     double *Q)  {
	    adapter::rysq::Shell A_(A), B_(B), C_(C), D_(D);
	    //rysq::Quartet<rysq::Shell> quartet(A_, B_, C_, D_);
	    boost::array<Center,4> centers = {{
		    A.center(), B.center(), C.center(), D.center() }};
	    rysq::Eri eri(rysq::Quartet<rysq::Shell>(A_, B_, C_, D_));
	    eri(centers, Q);
	}

	rysq_(const Shell::Data &A, const Shell::Data &B,
	      const Shell::Data &C, const Shell::Data &D)
	    : A_(A), B_(B), C_(C), D_(D),
	      eri_(rysq::Quartet<rysq::Shell>(A_, B_, C_, D_)) {}

	template<class A>
	void operator()(const std::vector<Center> &centers,
			const std::vector<boost::array<index,4> > &quartets,
			A &eri) {
	    BOOST_PROFILE_FUNCTION();
	    eri_(centers, quartets, eri);
	}

	void operator()(const std::vector<Center> &centers, size_t size,
			const boost::array<index,4> quartets[],
			double* const I[]) {
	    eri_(centers, size, quartets, I);
    }

		      
    private:
	adapter::rysq::Shell A_, B_, C_, D_;
	rysq::Eri eri_;
    };

    

    inline void evaluate(const Basis::Shell &A, const Basis::Shell &B,
    			 const Basis::Shell &C, const Basis::Shell &D,
			 double *Q) {
	rysq_::evaluate(A, B, C, D, Q);
    }

    inline void evaluate(const Basis::Shell &A, const Basis::Shell &B,
    			 const Basis::Shell &C, const Basis::Shell &D,
			 integral::array &Q) {
	size_t size = A.size()*B.size()*C.size()*D.size();
	if (Q.size() < size) {
	    std::cout << Q.size() << " " << size << std::endl;
	    throw std::length_error("insufficient size");
	}
	evaluate(A, B, C, D, Q.data());
    }

}

#endif // INTEGRAL_EVALUATE_HPP
