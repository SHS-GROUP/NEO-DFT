#ifndef INTEGRAL_HPP
#define INTEGRAL_HPP

#include "integral/forward.hpp"
#include "basis/basis.hpp"
#include "adapter/rysq.hpp"
#include "foreach.hpp"

#include "core/ublas.hpp"
#include <boost/array.hpp>

#include "core/integral/aux.hpp"
#include "core/integral/evaluate.hpp"
#include "core/integral/generator.hpp"

struct integral::Integral {

    typedef integral::index index;

    typedef integral::array Array;
    typedef integral::order Order;
    typedef integral::screening Screening;

    struct Gpu;

    Integral(const Basis::Shell &P, const Basis::Shell &Q,
	     const Basis::Shell &R, const Basis::Shell &S)
	: rysq_(P,Q,R,S) {}
    Integral(const Basis::Block &P, const Basis::Block &Q,
	     const Basis::Block &R, const Basis::Block &S)
	: rysq_(P.shell(), Q.shell(), R.shell(), S.shell()) {}

    // void operator()(const std::vector<Center> &centers,
    // 		    const std::vector<boost::array<index,4> > &quartets,
    // 		    const std::vector<double*> &I) {
    // 	rysq_(centers, quartets, I);
    // }
    // void operator()(const std::vector<Center> &centers,
    // 		    const std::vector<boost::array<index,4> > &quartets,
    // 		    double *I) {
    // 	rysq_(centers, quartets, I);
    // }
    void operator()(const std::vector<Center> &centers, size_t size,
		    const boost::array<index,4> quartets[],
		    double* const I[]) {
	rysq_(centers, size, quartets, I);
    }
private:
    integral::rysq_ rysq_;

public:
    struct Function {
	Function(const Basis &basis, const Screening &screening)
	    : basis_(basis), screening_(screening) {}

	template<typename T0, typename T1= void>
	struct result;

	template<typename T0>
	struct result<T0> { typedef integral::generator<1> type; };

	template<typename T0, typename T1>
	struct result { typedef integral::generator<2> type; };


	integral::generator<1> 
	operator()(size_t i) const {
	    return integral::make_generator(basis_, screening_, i);
	}

	integral::generator<2> 
	operator()(size_t i, size_t j) const {
	    return integral::make_generator(basis_, screening_, i, j);
	}

    private: 
	const Basis &basis_;
	const Screening &screening_;
    };

    static Function function(const Basis &basis, const Screening &screening) {
	return Function(basis, screening);
    }

};

typedef integral::Integral Integral;

#endif // INTEGRAL_HPP
