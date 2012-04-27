#ifndef INTEGRAL_ERI_RYSQ_HPP
#define INTEGRAL_ERI_RYSQ_HPP

#include "basis/basis.hpp"
#include "integral/quartet.hpp"

#include <boost/array.hpp>
#include <rysq.hpp>
#include "adapter/rysq.hpp"


namespace integral {
namespace eri {


    struct Rysq {

	Rysq(Quartet<const Basis::Shell::Data&> quartet) :
	    A_(quartet.p), B_(quartet.q),
	    C_(quartet.r), D_(quartet.s),
	    eri_(rysq::Quartet<rysq::Shell>(A_, B_, C_, D_)) {}

	template<class G>
	void operator()(const std::vector<Basis::Center> &centers,
			const std::vector< boost::array<int,4> > &quartets,
			G &g, double cutoff) {
	    //BOOST_PROFILE_FUNCTION();
	    // std::cout <<  "integrals: " << quartets.size() << std::endl;
	    eri_(centers, quartets, g, cutoff);
	}
		      
    private:
	adapter::rysq::Shell A_, B_, C_, D_;
	rysq::Eri eri_;
    };

}
}


#endif // INTEGRAL_ERI_RYSQ_HPP
