#ifndef DFT_GRID_HPP
#define DFT_GRID_HPP

#include "dft/dft.hpp"
#include "dft/grid/vectors.hpp"
#include "dft/functional.hpp"
#include "blas.hpp"
#include "dft/matrix.hpp"

#include "core/wavefunction.hpp"
#include "array_ref.hpp"

#define BOOST_PROFILE_ENABLE
#include "boost/utility/profiler.hpp"

#include <string>
#include <iostream>

#include <boost/typeof/typeof.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>



namespace dft {
namespace grid {

    struct Point {
	typedef boost::array<double,3> Center;

	template<class A>
	explicit Point(const A &dr) {
	    for (int i = 0; i < 3; ++i) {
		dr_[i] = dr[i];
	    }
	}
	Center operator()(const Basis::Shell &shell) const {
	    Center r = dr_;
	    for (int i = 0; i < 3; ++i) {
		r[i] += shell.center()[i];
	    }
	    return r;
	}
    private:
	Center dr_;
    };


}
}

namespace dft {

    struct Grid {

	typedef grid::Vectors Vectors;
	typedef grid::Point Point;
	typedef Functional::Density Density;
	typedef Functional::DE DE;

	explicit Grid(const Parameters &parameters = Parameters())
	    : p_(parameters) {}

	XC evaluate(const Wavefunction &W,
		    const_array_ref<double,3> dr,
		    const_array_ref<double> w,
		    const Functional &f) {

	    BOOST_PROFILE_FUNCTION();
	    
	    BOOST_AUTO(const &basis, W.basis());

	    vectors.resize_ao(basis.size(), dr.size());
	    vectors.mo().resize(W.occupied().size(), dr.size());
	    evaluate(W, dr, w);

	    if (f.gradient()) {
		vectors.resize_dq(basis.size(), dr.size());
		vectors.tmp().resize(basis.size(), dr.size());
		evaluate(W, dr);
	    }
	    
	    XC xc = f(density, de);

	    return xc;
	}

	template<class F_, class Side>
	void fock(boost::numeric::ublas::triangular_adaptor<F_,Side> F) {

	    BOOST_PROFILE_FUNCTION();

	    namespace blas = boost::numeric::bindings::blas;
	    namespace ublas = boost::numeric::ublas;

	    const size_t N = density.size();
	    int gradient = (vectors.dx().size1() ||
			    vectors.dy().size1() ||
			    vectors.dz().size1());

    	    ublas::range r1(0, vectors.ao().size1());
	    ublas::range r2(0, N);

    	    if (!gradient) {
    		for (size_t i = 0; i < N; ++i) {
    		    vectors.ao(i) *= sqrt(density.w[i]*de[0].dr[i]);
    		}
    		blas::syrk(1, ublas::project(vectors.ao(), r1, r2), 1, F);
    	    }
	    else {
		BOOST_AUTO(&tmp, vectors.tmp());
		BOOST_AUTO(X, ublas::project(tmp, r1, r2));
    		for (size_t i = 0; i < N; ++i) {    
		    double w = density.w[i];
		    ublas::column(X, i) =
			w*(de[0].dr[i]*vectors.ao(i)/2 +
			   de[0].dx[i]*vectors.dx(i) +
			   de[0].dy[i]*vectors.dy(i) +
			   de[0].dz[i]*vectors.dz(i));
    		}
		// blas::syr2k(1, ublas::project(vectors.ao(), r1, r2), X,
		// 	     1, F);
		block::outer_prod<16>
		    (1.0, ublas::project(vectors.ao(), r1, r2), X,
		     1.0, F, p_.ccut);
    	    }

	    if (gradient == 2) {
    		// blas::syrk(gga, project(vectors.dx(), r1, r2), F);
    		// blas::syrk(gga, project(vectors.dy(), r1, r2), F);
    		// blas::syrk(gga, project(vectors.dz(), r1, r2), F);
    	    }

    	}

    private:
	Parameters p_;
	Vectors vectors;
	Density density;
	DE de;

	void evaluate(const Wavefunction &W,
		      const_array_ref<double,3> dr,
		      const_array_ref<double> w);

	void evaluate(const Wavefunction &W,
		      const_array_ref<double,3> dr);

    };


}

#endif // DFT_GRID_HPP
