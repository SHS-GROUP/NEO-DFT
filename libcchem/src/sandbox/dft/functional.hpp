#ifndef DFT_FUNCTIONAL_HPP
#define DFT_FUNCTIONAL_HPP

#include "dft/xc/functional.hpp"

#include <memory>
#include <vector>

#include <boost/noncopyable.hpp>
#include <boost/typeof/typeof.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <boost/numeric/bindings/blas.hpp>
#include <boost/numeric/bindings/ublas.hpp>

namespace dft {    

    struct Functional : boost::noncopyable {
	typedef xc::functional::Density Density;
	typedef xc::functional::XC XC;
	typedef xc::functional::DE DE;
	explicit Functional(const std::string &name)
	    : impl_(xc::functional::new_(name)) {}
	int gradient() const { return impl_->gradient; }
	XC operator()(const Density &density, DE &de) const {
	    XC xc = { 0 };
	    apply(*impl_, density, de, xc);
	    return xc;
	}
    private:
	std::auto_ptr<xc::functional> impl_;
	template<class F, class R, class D, class X>
	static void apply(F &f, const R &density, D &de, X &xc) {
	    f(density, de, xc);
	}
    };


}


#endif // DFT_FUNCTIONAL_HPP
