#ifndef DFT_XC_GAMESS_HPP
#define DFT_XC_GAMESS_HPP

#include "dft/xc/kernel.hpp"

struct gamess {

    template<class F>
    static void exchange(const F &f, double scale,
			 const Point::Density &density,
			 Point::XC &xc) {
	double ignore = 0;
	f(scale, true, false,
	  density.a, density.b,
	  density.aa, density.bb,
	  density.dx, density.dy, density.dz,
	  0, 0, 0,
	  density.w,
	  xc.sigk, ignore,
	  xc.dk, ignore,
	  xc.Xa, xc.Xg,
	  xc.V, xc.dx, xc.dy, xc.dz,
	  ignore, ignore, ignore, ignore);
    }

    template<class F>
    static void correlation(const F &f, double scale,
			 const Point::Density &density,
			 Point::XC &xc) {
	double ignore = 0;
	f(scale, true, false, 0,
	  density.a, density.b,
	  density.aa, density.bb, density.ab,
	  density.dx, density.dy, density.dz,
	  0, 0, 0,
	  density.w,
	  xc.Ec,
	  xc.V, xc.dx, xc.dy, xc.dz,
	  ignore, ignore, ignore, ignore);
    }

};

#endif // DFT_XC_GAMESS_HPP
