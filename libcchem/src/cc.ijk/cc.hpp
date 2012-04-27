#ifndef CC_CC_HPP
#define CC_CC_HPP

#include "array/array.hpp"
#include <cstdlib>
#include <map>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/ref.hpp>

namespace cc {

    struct Triples {
	typedef boost::numeric::ublas::vector<double> vector;


	struct Correction {
	    double ots, otd, ets, etd;
	    Correction() : ots(0), otd(0), ets(0), etd(0) {}
	    Correction& operator+=(const Correction &c) {
		ots += c.ots;
		otd += c.otd;
		ets += c.ets;
		etd += c.etd;
		return *this;
	    }
	};

	typedef std::map<
	    std::string,
	    boost::reference_wrapper<const Array<double> >
	    > Data;


	Correction operator()(size_t no, size_t nv,
			      const Data &data,
			      const vector &eh, const vector &ep,
			      double* t3);

	// Correction operator()(size_t no, size_t nv,
	// 		      const Data &data,
	// 		      const vector &eh, const vector &ep,
	// 		      int i, int j, int k,
	// 		      double* u1, double* t3);
	
    };

}

#endif // CC_CC_HPP
