#ifndef INTEGRAL_SCREENING_HPP
#define INTEGRAL_SCREENING_HPP

#include "basis/basis.hpp"

#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include "foreach.hpp"

namespace integral {

    struct Screening {

	typedef boost::numeric::ublas::matrix<
	    double,
	    boost::numeric::ublas::column_major
	    > Matrix;

	const double value;

	Screening(size_t n) :
	    value(0), K_(n,n) 
	{
	    std::fill(K_.data().begin(), K_.data().end(), 1);
	}

	template<class E>
	Screening(const E &e, double value) :
	    value(value), K_(e) {}

	template<class I>
	void operator()(const Basis::Block &P, const Basis::Block &Q,
			const Basis::Block &R, const Basis::Block &S,
			std::vector<I> &quartets) const {
	    evaluate(P, Q, R, S, quartets, *this);
	}
	template<class Q>
	bool operator()(const Q &q) const {
	    return !(K_(q[0],q[1])*K_(q[2],q[3]) < value);
	}

	template<class I, class T>
	static void evaluate(const Basis::Block &P, const Basis::Block &Q,
			     const Basis::Block &R, const Basis::Block &S,
			     std::vector<I> &quartets, T &test) {
	    typedef Basis::Shell Shell;
	    foreach (const Shell &s, S) {
		foreach (const Shell &r, R) {
		    foreach (const Shell &q, Q) {
			foreach (const Shell &p, P) {
			    I quartet =
				{{ p.index(), q.index(), r.index(), s.index() }};
			    if (!test(quartet)) continue;
			    quartets.push_back(quartet);
			}
		    }
		}
	    }
	}

    private:
	Matrix K_;
    };

}

#endif // INTEGRAL_SCREENING_HPP
