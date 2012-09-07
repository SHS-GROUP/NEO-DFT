#include "dft/xc/kernel.hpp"
#include "dft/xc/gamess.hpp"

struct b88 : needs_gradient<1> {

    const double scale_;
    explicit b88(double scale) : scale_(scale) {}

    void operator()(const Point::Density &density, Point::XC &xc) {
	gamess::exchange(&b88xf_, scale_, density, xc);
    }

    static
    int b88xf_(const double scale,
	       const bool rhfgvb, const bool urohf,
	       const double roa, const double rob,
	       const double grdaa, const double grdbb,
	       const double gradxa, const double gradya, const double gradza,
	       const double gradxb, const double gradyb, const double gradzb,
	       const double ftotwt,
	       double &sigka, double &sigkb,
	       double &dkdxa, double &dkdxb,
	       double &xalf, double &b88x,
	       double &duma, double &dumxa, double &dumya, double &dumza,
	       double &dumb, double &dumxb, double &dumyb, double &dumzb) {

	static const double c_b2 = .33333333333333331;
	static const double c_b3 = .5;
	static const double c_b4 = .66666666666666663;

	/* System generated locals */
	double d__1, d__2;

	/* Local variables */
	double e, g, h__, q, r__, s, t, u, w, x, x1, ea, eb, ga, ha, gb,
	    hb, qa, ra, sa, ta, ua, qb, wa, xa, xb, sb, tb, rb, ub, wb,
	    grd, xcl, rho, face, facp, xcla, xclb, b88xa, rho13, rhoe,
	    b88xb, xcnl, xalfa, xalfb, rho13a, rhoea, rho13b, pduma, pdumb,
	    rhoeb, xcnla, xcnlb, beckeb;
	double pdumxa, pdumya, pdumza, pdumxb, pdumyb, pdumzb;

	/* *********************************************************************** */
	/*     BECKE 88 EXCHANGE FUNCTIONAL */

	/*     A. D. BECKE */
	/*     PHYS. REV. A 38,3098 (1988) */

	/*     ----- */
	/*     INPUT */
	/*     ----- */

	/*     ROA,B   .... SPINDENSITIES */
	/*     GRAD**  .... DENSITY GRADIENTS */
	/*     FTOTWT  .... VOLUME OF THIS GRID */

	/*     ------ */
	/*     OUTPUT */
	/*     ------ */

	/*     XALF    ....  EXCHANGE ENERGY OF LDA PART */
	/*     B88X    ....  EXCHANGE ENERGY OF GGA PART */
	/*     DUM***  ....  FIRST DERIVATIVES OF THE EXCHANGE FUNCTINALS */
	/*     SIGKA,B ....  EXC=-1/2*INT(SIGK*RHO^(4/3)) */
	/*     DKDXA,B ....  DSIGK/DX */
	/* *********************************************************************** */
	/* *********************************************************************** */
	d__1 = 3. / PI;
	x1 = pow(d__1, c_b2);
	beckeb = .0042;
	face = x1*-1.5*pow(c_b3, c_b4);
	facp = 1.3333333333333333;
	/* *********************************************************************** */
	/*      CLOSED SHELL CALCULATION */
	/* *********************************************************************** */
	if (rhfgvb) {
	    rho = roa;
	    grd = grdaa;
	    rho13 = pow(rho, c_b2);
	    /* Computing 4th power */
	    d__1 = rho13, d__1 *= d__1;
	    rhoe = d__1*d__1;
	    x = sqrt(grd) / rhoe;
	    q = sqrt(x*x + 1.);
	    s = log(x + q);
	    t = beckeb*6.*x*s + 1.;
	    g = -beckeb*x*x / t;
	    /* Computing 2nd power */
	    d__1 = beckeb*x;
	    /* Computing 2nd power */
	    d__2 = t;
	    h__ =
		(beckeb*-2.*x +
		 d__1*d__1*6.*(-s + x / q)) / (d__2*d__2);
	    r__ = facp*rho13;
	    e = face;
	    xcl = e*rhoe;
	    xcnl = g*rhoe;
	    xalfa = scale*2.*xcl*ftotwt;
	    b88xa = scale*2.*xcnl*ftotwt;
	    u = scale*r__*(e + g - x*h__);
	    pduma = u;
	    w = scale*.5*h__ / sqrt(grd);
	    pdumxa = w*2.*gradxa;
	    pdumya = w*2.*gradya;
	    pdumza = w*2.*gradza;
	    /*     ----- FOR OP CORRELATION ----- */
	    sigka = -(face + g)*2.;
	    dkdxa = -h__*2.;
	    /*     ------------------------------ */
	    // if (nlrc_1.lcflag) {
	    // 	lrcfun_(urohf, roa, rob, grdaa, grdbb, gradxa, gradya,
	    // 		gradza, gradxb, gradyb, gradzb, ftotwt, sigka,
	    // 		sigkb, dkdxa, dkdxb, xalfa, xalfb, b88xa,
	    // 		b88xb, pduma, pdumxa, pdumya, pdumza, pdumb,
	    // 		pdumxb, pdumyb, pdumzb, face);
	    // }
	    xalf += xalfa;
	    b88x += b88xa;
	    duma += pduma;
	    dumxa += pdumxa;
	    dumya += pdumya;
	    dumza += pdumza;
	    /* *********************************************************************** */
	    /*      OPEN SHELL CALCULATION */
	    /* *********************************************************************** */
	}
	else if (urohf) {
	    /*     ----- ALPHA CONTRIBTUTION ---- */
	    rho13a = pow(roa, c_b2);
	    /* Computing 4th power */
	    d__1 = rho13a, d__1 *= d__1;
	    rhoea = d__1*d__1;
	    xa = sqrt(grdaa) / rhoea;
	    /* Computing 2nd power */
	    d__1 = xa;
	    qa = sqrt(d__1*d__1 + 1.);
	    sa = log(xa + qa);
	    ta = beckeb*6.*xa*sa + 1.;
	    ga = -beckeb*xa*xa / ta;
	    /* Computing 2nd power */
	    d__1 = beckeb*xa;
	    /* Computing 2nd power */
	    d__2 = ta;
	    ha = (beckeb*-2.*xa +
		  d__1*d__1*6.*(-sa + xa / qa)) / (d__2*d__2);
	    ra = facp*rho13a;
	    ea = face;
	    ua = scale*ra*(ea + ga - xa*ha);
	    pduma = ua;
	    wa = scale*.5*ha / sqrt(grdaa);
	    pdumxa = wa*2.*gradxa;
	    pdumya = wa*2.*gradya;
	    pdumza = wa*2.*gradza;
	    /*     ----- BETA CONTRIBTUTION ---- */
	    rho13b = pow(rob, c_b2);
	    /* Computing 4th power */
	    d__1 = rho13b, d__1 *= d__1;
	    rhoeb = d__1*d__1;
	    eb = face;
	    if (rob < 1e-15) {
		gb = 0.;
		hb = 0.;
		pdumb = 0.;
		pdumxb = 0.;
		pdumyb = 0.;
		pdumzb = 0.;
		goto L9999;
	    }
	    xb = sqrt(grdbb) / rhoeb;
	    /* Computing 2nd power */
	    d__1 = xb;
	    qb = sqrt(d__1*d__1 + 1.);
	    sb = log(xb + qb);
	    tb = beckeb*6.*xb*sb + 1.;
	    gb = -beckeb*xb*xb / tb;
	    /* Computing 2nd power */
	    d__1 = beckeb*xb;
	    /* Computing 2nd power */
	    d__2 = tb;
	    hb = (beckeb*-2.*xb +
		  d__1*d__1*6.*(-sb + xb / qb)) / (d__2*d__2);
	    rb = facp*rho13b;
	    ub = scale*rb*(eb + gb - xb*hb);
	    pdumb = ub;
	    wb = scale*.5*hb / sqrt(grdbb);
	    pdumxb = wb*2.*gradxb;
	    pdumyb = wb*2.*gradyb;
	    pdumzb = wb*2.*gradzb;
	L9999:
	    xcla = ea*rhoea;
	    xclb = eb*rhoeb;
	    xcnla = ga*rhoea;
	    xcnlb = gb*rhoeb;
	    xalfa = scale*xcla*ftotwt;
	    xalfb = scale*xclb*ftotwt;
	    b88xa = scale*xcnla*ftotwt;
	    b88xb = scale*xcnlb*ftotwt;
	    /*     ----- FOR OP CORRELATION ----- */
	    sigka = -(face + ga)*2.;
	    sigkb = -(face + gb)*2.;
	    dkdxa = -ha*2.;
	    dkdxb = -hb*2.;
	    /*     ------------------------------- */
	    // if (nlrc_1.lcflag) {
	    // 	lrcfun_(urohf, roa, rob, grdaa, grdbb, gradxa, gradya,
	    // 		gradza, gradxb, gradyb, gradzb, ftotwt, sigka,
	    // 		sigkb, dkdxa, dkdxb, xalfa, xalfb, b88xa,
	    // 		&b88xb, pduma, pdumxa, pdumya, pdumza, pdumb,
	    // 		&pdumxb, pdumyb, pdumzb, face);
	    // }
	    xalf = xalf + xalfa + xalfb;
	    b88x = b88x + b88xa + b88xb;
	    duma += pduma;
	    dumxa += pdumxa;
	    dumya += pdumya;
	    dumza += pdumza;
	    dumb += pdumb;
	    dumxb += pdumxb;
	    dumyb += pdumyb;
	    dumzb += pdumzb;
	}
	return 0;
    }

};
