#ifndef CC_CC_HPP
#define CC_CC_HPP

#include "core/wavefunction.hpp"
#include "array/array.hpp"

#include <cstdlib>
#include <map>
#include <stdexcept>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/ref.hpp>

namespace cc {

    template<class P>
    struct Map {
	struct no_such_key : std::runtime_error {
	    no_such_key(const char *what) : std::runtime_error(what) {}
	};
	P& operator()(const char *key) {
	    // std::cout << "[" << key << "]" << std::endl;
	    return data_[key];
	}
	const P& operator()(const char *key) const {
	    // std::cout << "(" << key << ")" << std::endl;
	    if (!has(key)) throw no_such_key(key);
	    return data_.find(key)->second;
	}
	bool has(const char *key) const {
	    return (data_.find(key) != data_.end());
	}
	const std::map<std::string, P>& data() const {
	    return data_;
	}
    private:
	std::map<std::string, P> data_;
    };

    template<class P>
    std::ostream& operator<<(std::ostream &os, const Map<P> &m) {
	BOOST_AUTO(it, m.data().begin());
	while (it != m.data().end()) {
	    os << "[" << it->first << "]";
	    ++it;
	}
	return os;
    }

    typedef boost::numeric::ublas::vector<double> Vector;
    typedef boost::numeric::ublas::column_major Layout;
    typedef boost::numeric::ublas::matrix<double, Layout> Matrix;

    typedef ::Array<double> Array;

}

namespace cc {
namespace sd {

    void vvvv(const Wavefunction &wf, Map<const Array*> t,
	      Array &u2, Array &c);

}
}

namespace cc {

    struct Triples {

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
	    boost::reference_wrapper<const Array>
	    > Data;


	Correction operator()(size_t no, size_t nv,
			      Map<const Array*> t,
			      Map<const Array*> V,
			      const Map<Vector> &e);

	// Correction operator()(size_t no, size_t nv,
	// 		      const Data &data,
	// 		      const vector &eh, const vector &ep,
	// 		      int i, int j, int k,
	// 		      double* u1, double* t3);
	
    };

    typedef Triples::Correction Correction;

}


#endif // CC_CC_HPP
