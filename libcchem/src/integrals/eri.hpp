#ifndef INTEGRALS_ERI_HPP
#define INTEGRALS_ERI_HPP

#include "basis/basis.hpp"
#include "integrals/quartet.hpp"
#include "integrals/screening.hpp"
#include "integrals/rysq.hpp"

#include "boost/utility/profiler.hpp"

#include <vector>
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "foreach.hpp"

namespace integrals {
namespace detail {

    inline Basis::Block block(const Basis::Block &b) {
    	return b;
    }

    inline  Basis::Block block(const Basis::Shell &s) {
	return Basis::Block(&s);
    }

    inline const Basis::Shell::Data& shell(const Basis::Shell &s) {
	return *s.data();
    }

    inline const Basis::Shell::Data& shell(const Basis::Block &b) {
    	return b.shell();
    }


}
}

namespace integrals {

    struct Eri {

	struct Cuda;

	typedef Basis::Shell Shell;
	typedef Basis::Center Center;
	typedef std::vector< boost::array<int,4> > Quartets;

	Eri(const Basis &basis, const Screening *screen)
	    : centers_(basis.centers()), screen_(screen)
	{
	    max_ = basis.data().size();
	}

	template<class P, class Q, class R, class S>
	void screen(const P &p, const Q &q, const R &r, const S &s,
		    Quartets &quartets) const {
	    //BOOST_PROFILE_LINE;
	    (*screen_)(detail::block(p), detail::block(q),
		       detail::block(r), detail::block(s),
		       quartets);
	}

	template<class P, class Q, class R, class S>
	void operator()(const P &p, const Q &q, const R &r, const S &s,
			const Quartets &quartets, double *data) {
	    //BOOST_PROFILE_LINE;
	    Quartet<const Shell::Data&> quartet =
		{ detail::shell(p), detail::shell(q),
		  detail::shell(r), detail::shell(s) };
	    std::fill_n(data, quartets.size()*quartet.size(), 0);
	    evaluate(quartet, centers_, quartets, data, screen_->value());
	}

	template<class P, class Q, class R, class S>
	void operator()(const P &p, const Q &q, const R &r, const S &s,
			double *data = NULL) {
	    //BOOST_PROFILE_LINE;

	    Quartet<const Shell::Data&> quartet =
		{ detail::shell(p), detail::shell(q),
		  detail::shell(r), detail::shell(s) };

	    {
		quartets_.clear();
		(*screen_)(detail::block(p), detail::block(q),
			   detail::block(r), detail::block(s),
			   quartets_);
	    }

	    if (!data) {
		data_.reserve(quartets_.size()*quartet.size());
		data = &data_[0];
	    }
	    std::fill_n(data, quartets_.size()*quartet.size(), 0);

	    evaluate(quartet, centers_, quartets_, data, screen_->value());
	}

	const Quartets& quartets() const {
	    return quartets_;
	}

	const double* data() const { return &data_[0]; }

    private:
	std::vector<Basis::Center> centers_;
	Quartets quartets_;
	std::vector<double> data_;
	const Screening *screen_;

	typedef boost::array<Shell::Data::key_type,4> key_type;
	boost::ptr_map<key_type, integrals::rysq::Eri> cache_;
	size_t max_;

	void evaluate(Quartet<const Shell::Data&> quartet,
		      const std::vector<Center> &centers,
		      const Quartets &quartets,
		      double *G, double cutoff) {
	    //BOOST_PROFILE_LINE;
	    key_type key = {{ quartet.p.key(), quartet.q.key(),
			      quartet.r.key(), quartet.s.key() }};
	    if (cache_.size() > max_) {
		cache_.clear();
	    }
	    if (!cache_.count(key)) {
		BOOST_PROFILE_LINE;
		// std::cout << "insert: "
		// 	  << key[0] << " " << key[1] << " "
		// 	  << key[2] << " " << key[3] 
		// 	  << std::endl;
		cache_.insert(key, new integrals::rysq::Eri(quartet));
	    }
	    integrals::rysq::Eri &rysq = cache_.at(key);
	    rysq(centers, quartets, G, cutoff);
	}

    private:

    };

#ifdef CCHEM_INTEGRALS_ERI_CUDA
    struct Eri::Cuda {
	Cuda(const Basis &basis, const Screening *screen)
	    : rysq_(basis), screen_(screen) {}
	template<class P, class Q, class R, class S>
	void operator()(const P &p, const Q &q, const R &r, const S &s,
			const Quartets &quartets, double *data) {
	    cuda::stream stream(0);
	    //BOOST_PROFILE_LINE;
	    // std::fill_n(data, quartets.size()*quartet.size(), 0);
	    rysq_(detail::shell(p),
	    	  detail::shell(q),
	    	  detail::shell(r),
	    	  detail::shell(s),
	    	  quartets, data,
	    	  screen_->value(),
	    	  stream);
	    stream.synchronize();
	}
    private:
	rysq::Eri::Cuda rysq_;
	const Screening *screen_;
    };
#endif // CCHEM_INTEGRALS_ERI_CUDA

}

#endif // INTEGRALS_ERI_HPP
