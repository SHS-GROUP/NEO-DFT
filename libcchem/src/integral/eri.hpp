#ifndef INTEGRAL_ERI_HPP
#define INTEGRAL_ERI_HPP

#include "basis/basis.hpp"
#include "integral/quartet.hpp"
#include "integral/screening.hpp"
#include "integral/eri/rysq.hpp"

#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include "foreach.hpp"

namespace integral {
namespace detail {


    static Basis::Block block(const Basis::Block &b) {
	return b;
    }

    static Basis::Block block(const Basis::Shell &s) {
	return Basis::Block(&s);
    }

    static const Basis::Shell::Data& shell(const Basis::Shell &s) {
	return *s.data();
    }

    static const Basis::Shell::Data& shell(const Basis::Block &b) {
	return b.shell();
    }


}
}

namespace integral {


    struct Eri {

	typedef Basis::Shell Shell;
	typedef Basis::Center Center;
	typedef boost::array<int,4> int4;

	Eri(const Basis &basis, const Screening *screen)
	    : centers_(basis.centers()), screen_(screen) {}

	template<class P, class Q, class R, class S>
	void operator()(const P &p, const Q &q, const R &r, const S &s) {
	    Quartet<const Shell::Data&> quartet =
		{ detail::shell(p), detail::shell(q),
		  detail::shell(r), detail::shell(s) };

	    quartets_.clear();
	    (*screen_)(detail::block(p), detail::block(q),
		       detail::block(r), detail::block(s),
		       quartets_);

	    data_.reserve(quartets_.size()*quartet.size());
	    data_.clear();
	    data_.resize(quartets_.size()*quartet.size(), 0);
	    evaluate(quartet, centers_, quartets_, &data_[0], screen_->value);
	}

	void evaluate(Quartet<const Shell::Data&> quartet,
		      const std::vector<Center> &centers,
		      const std::vector<int4> &quartets,
		      double *G, double cutoff) {
	    integral::eri::Rysq rysq(quartet);
	    rysq(centers, quartets, G, cutoff);
	}

	const std::vector<int4>& quartets() const {
	    return quartets_;
	}

	const double* data() const { return &data_[0]; }

    private:
	std::vector<Basis::Center> centers_;
	std::vector<int4> quartets_;
	std::vector<double> data_;
	const Screening *screen_;

    private:

    };

}

#endif // INTEGRAL_ERI_HPP
