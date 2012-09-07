#ifndef DFT_XC_FUNCTIONAL_HPP
#define DFT_XC_FUNCTIONAL_HPP

#include <string>
#include <cmath>
#include <vector>
#include <iostream>

#include <boost/array.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/functional.hpp>

#include "dft/types.hpp"

namespace dft {
namespace xc {


#include "dft/xc/kernel.hpp"


    template<class F>
    struct kernel;

    struct functional {
	class Density;
	class DE;
	typedef dft::XC XC;

	const int gradient;
	functional(int gradient) : gradient(gradient) {};
	virtual ~functional() {}

	virtual
	void operator()(const Density &density, DE &de, XC &xc) const = 0;

	static functional* new_(const std::string &name);
	template<class F>
	static functional* new_(const F &f) {
	    return new kernel<F>(f);
	}
    };


    struct functional::Density {
	typedef std::vector<double> vector_type;

	std::vector<int> index;
	vector_type w, a, b;
	vector_type dx, dy, dz;
	vector_type aa, bb, ab;

	size_t size() const { return index.size(); }
	void push_back(int i, double w, const double (&rho)[2], int npol) {
	    this->index.push_back(i);
	    this->w.push_back(w);
	    this->a.push_back(rho[0]);
	    if (npol == 2) this->b.push_back(rho[1]);
	}

	Point::Density operator[](int i) const {
	    Point::Density p = { w[i],
				 a[i], a[i], aa[i], aa[i], aa[i],
				 dx[i], dy[i], dz[i] };
	    if (int(b.size()) > i) {
		p.b = b[i];
		p.bb = bb[i];
		p.ab = ab[i];
	    }
	    return p;
	}
    };

    struct functional::DE {
	struct tuple {
	    std::vector<double> dr, dx, dy, dz;
	    void resize(size_t size) {
		dr.resize(size);
		dx.resize(size);
		dy.resize(size);
		dz.resize(size);
	    }
	    struct reference {
		typedef double& type;
		type dr, dx, dy, dz;
		void operator=(const Point::XC &p) {
		    dr = p.V;
		    dx = p.dx;
		    dy = p.dy;
		    dz = p.dz;
		}
	    };
	    reference operator[](int i) {
		reference r = { dr[i], dx[i], dy[i], dz[i] };
		return r;
	    }
	};
	tuple& operator[](bool i) { return data_.at(i); }
    private:
	boost::array<tuple,2> data_;
    };


    template<class F>
    struct kernel : functional {
	F f;
	explicit kernel(const F &f, int dr = 0)
	    : functional(max_gradient(f, dr).value), f(f) {}
	void operator()(const Density &density, DE &de, XC &xc) const {
	    namespace fusion = boost::fusion;
	    de[0].resize(density.size());

#if (defined(__GNUC__) && ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3)))
#warning "GNUC < 4.3 with -fopenmp causes ICE, needs fix"
#else
#pragma omp parallel
#endif
	    {
		XC xc_ = { 0 };
		typedef fusion::vector<const Point::Density&, Point::XC&> A;
#pragma omp for schedule(dynamic,16)
		for (int i = 0; i < int(density.size()); ++i) {
		    const Point::Density& d = density[i];
		    Point::XC p = { 0 };
		    boost::fusion::for_each(f, apply<A>(A(d, p)));
		    de[0][i] = p;
		    xc_.Xa += p.Xa;
		    xc_.Xg += p.Xg;
		    xc_.Ec += p.Ec;
		    xc_.rho += d.w*(d.a+d.b);
		    xc_.aa += d.w*d.aa;
		    xc_.bb += d.w*d.bb;
		    xc_.ab += d.w*d.ab;
		}
		#pragma omp critical
		xc += xc_;
	    }
	}
    private:
	struct max_gradient {
	    mutable int value;
	    max_gradient(const F &f, int value) {
		this->value = value;
		boost::fusion::for_each(f, *this);
	    }
	    BOOST_MPL_HAS_XXX_TRAIT_DEF(gradient);
	    template<class F_>
	    typename boost::enable_if<has_gradient<F_> >::type
	    operator()(const F_ &f) const {
		static const int g = F_::gradient::value;
	    	if (g > value) value = g;
	    }
	    template<class F_>
	    typename boost::disable_if<has_gradient<F_> >::type
	    operator()(const F_ &f) const {}
	};
	template<class A>
	struct apply {
	    A args;
	    explicit apply(const A &args) : args(args) {}
	    template<class F_>
	    void operator()(const F_ &f) const {
		boost::fusion::invoke_procedure(f, args);
	    }
	};
    };
    
}
}

#endif // DFT_XC_FUNCTIONAL_HPP

