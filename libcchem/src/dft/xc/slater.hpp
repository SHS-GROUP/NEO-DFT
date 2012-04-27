#include "dft/xc/kernel.hpp"

struct slater {

    const double scale;
    explicit slater(double scale) : scale(scale) {}

    void operator()(const Point::Density &density, Point::XC &xc) const {
	double ignore;
    	lsdslt_(true, false, density.a, density.b, density.w,
    		xc.V, ignore, xc.Xa, xc.sigk, ignore);
    }

    int lsdslt_(bool rhfgvb, bool urohf,
		const double roa, const double rob, const double ftotwt,
		double &vxca, double &vxcb, double &xalpha,
		double &sigka, double &sigkb) const {

	static double c_b2 = .33333333333333331;

	/* System generated locals */
	double d__1;

	/* Local variables */
	double face, facp, rho13, rhoe, grdaa, grdbb,
	    dkdxa, dkdxb, xalfa, xalfb, rho13a, rho13b, rhoea, pduma,
	    pdumb, rhoeb, gradxa, gradya, gradza, gradxb, gradyb, gradzb;

	//double b88xa, b88xb, pdumxa, pdumya, pdumza, pdumxb, pdumyb, pdumzb;

	/* *********************************************************************** */
	/*     SLATER(DIRAC) EXCHANGE FUNCTIONAL */

	/*     XALF    ....  EXCHANGE ENERGY OF LDA PART */
	/*     B88X    ....  EXCHANGE ENERGY OF GGA PART */
	/*     DUM     ....  EXCHANGE GRADIENT  TO RHO */
	/*     DUMX    ....  EXCHANGE GRADIENT TO X */
	/*     DUMY    ....  EXCHANGE GRADIENT TO Y */
	/*     DUMZ    ....  EXCHANGE GRADIENT TO Z */
	/* *********************************************************************** */
	/* *********************************************************************** */
	grdaa = 0.;
	grdbb = 0.;
	gradxa = 0.;
	gradya = 0.;
	gradza = 0.;
	gradxb = 0.;
	gradyb = 0.;
	gradzb = 0.;
	dkdxa = 0.;
	dkdxb = 0.;
	if (rhfgvb) {
	    rho13 = pow(roa, c_b2);
	    d__1 = 6. / PI;
	    facp = -pow(d__1, c_b2);
	    d__1 = 6. / PI;
	    face = pow(d__1, c_b2)*-.75;
	    pduma = facp*rho13;
	    vxcb = vxcb;
	    /* Computing 4th power */
	    d__1 = rho13, d__1 *= d__1;
	    rhoe = d__1*d__1;
	    xalfa = ftotwt*rhoe*face*2.;
	    /*     ----- FOR OP CORRELATION ----- */
	    sigka = -face*2.;
	    /*     ------------------------------ */
	    // if (nlrc_1.lcflag) {
	    //     dkdxa = 0.;
	    //     lrcfun_(urohf, roa, rob, grdaa, grdbb, gradxa, gradya, 
	    //          gradza, gradxb, gradyb, gradzb, ftotwt, sigka, sigkb, 
	    //          dkdxa, dkdxb, xalfa, xalfb, b88xa, b88xb, pduma, 
	    //          pdumxa, pdumya, pdumza, pdumb, pdumxb, pdumyb, 
	    //          pdumzb, face);
	    // }
	    xalpha += scale*xalfa;
	    vxca += scale*pduma;
	}
	else if (urohf) {
	    rho13a = pow(roa, c_b2);
	    rho13b = pow(rob, c_b2);
	    d__1 = 6. / PI;
	    facp = -pow(d__1, c_b2);
	    d__1 = 6. / PI;
	    face = pow(d__1, c_b2)*-.75;
	    pduma = scale*facp*rho13a;
	    pdumb = scale*facp*rho13b;
	    /* Computing 4th power */
	    d__1 = rho13a, d__1 *= d__1;
	    rhoea = d__1*d__1;
	    /* Computing 4th power */
	    d__1 = rho13b, d__1 *= d__1;
	    rhoeb = d__1*d__1;
	    xalfa = scale*ftotwt*rhoea*face;
	    xalfb = scale*ftotwt*rhoeb*face;
	    /*     ----- FOR OP CORRELATION ----- */
	    sigka = -face*2.;
	    sigkb = -face*2.;
	    /*     ------------------------------ */
	    // if (nlrc_1.lcflag) {
	    // 	dkdxa = 0.;
	    // 	dkdxb = 0.;
	    /*     PERHAPS NEED NOT GRDAA,GRDBB=0.0D+00 */
	    // lrcfun_(urohf, roa, rob, grdaa, grdbb, gradxa, gradya, 
	    //      gradza, gradxb, gradyb, gradzb, ftotwt, sigka, sigkb, 
	    //      dkdxa, dkdxb, xalfa, xalfb, b88xa, b88xb, pduma, 
	    //      pdumxa, pdumya, pdumza, pdumb, pdumxb, pdumyb, 
	    //      pdumzb, face);
	    // }
	    xalpha = xalpha + xalfa + xalfb;
	    vxca += pduma;
	    vxcb += pdumb;
	}
	return 0;
    }

};
