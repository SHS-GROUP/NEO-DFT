#ifndef INTEGRALS_SCREENING_HPP
#define INTEGRALS_SCREENING_HPP

#include "basis/basis.hpp"

#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "foreach.hpp"

namespace integrals {

    struct Screening {

	typedef boost::numeric::ublas::matrix<
	    double,
	    boost::numeric::ublas::column_major
	    > Matrix;

    private:
	Matrix K_;
	double value_;

    public:

	Screening(size_t n) :
	    K_(boost::numeric::ublas::scalar_matrix<double>(n,n,0)) ,
	    value_(0) {}

	template<class E>
	Screening(const E &e, double value) :
	    K_(e), value_(value) {}

	Screening(const Basis &basis, double value);

	double value() const { return value_; }

	template<class I>
	void operator()(const Basis::Block &P, const Basis::Block &Q,
			const Basis::Block &R, const Basis::Block &S,
			std::vector<I> &quartets) const {
	    evaluate(P, Q, R, S, quartets, *this);
	}
	template<class Q>
	bool operator()(const Q &q) const {
	    // std::cout << "yyy "
	    // 	      << q[0] <<"," << q[1] <<","
	    // 	      << q[2] <<"," << q[3] <<" " 
	    // 	      << K_(q[0],q[1]) << ","
	    // 	      << K_(q[2],q[3]) << " = "
	    // 	      << (K_(q[0],q[1])*K_(q[2],q[3])) << " -> "
	    // 	      << !(K_(q[0],q[1])*K_(q[2],q[3]) < value_) <<  std::endl;
	    return !(K_(q[0],q[1])*K_(q[2],q[3]) < value_);
	}

	template<class V>
	void reorder(const V &index) {
	    namespace ublas = boost::numeric::ublas;
	    ublas::indirect_array<> p(index.size());
	    for (size_t i = 0; i < index.size(); ++i) {
		p[i] = *(index.begin()+i);
	    }
	    K_.assign(Matrix(ublas::project(K_, p, p)));
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
    };

}

#endif // INTEGRALS_SCREENING_HPP
