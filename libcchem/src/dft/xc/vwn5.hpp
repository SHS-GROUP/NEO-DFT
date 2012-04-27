#include "dft/xc/kernel.hpp"

struct vwn5 {

    const double scale_;
    explicit vwn5(double scale) : scale_(scale) {}

    void operator()(const Point::Density &density, Point::XC &xc) {
	double ignore = 0;
	vwn5cf_(scale_, true, false,
		density.a, density.b, density.w,
		xc.V, ignore, xc.Ec);
    }

    static
    int vwn5cf_(const double scale,
		const bool rhfgvb, const bool urohf,
		const double roa, const double rob, const double ftotwt, 
		double &vxca, double &vxcb, double &evwn) {

	static const double c_b2 = .16666666666666666;
	static const double c_b4 = 1.3333333333333333;
	static const double c_b6 = .33333333333333331;

	/* System generated locals */
	double d__1, d__2;

	/* Local variables */
	double q, x, t1, t2, t3, qa, qf, hx, qp, xx, t1a, t2a, t3a,
	    t1f, t2f, t3f, t1p, t2p, t3p, xx0, xba, xbf, dhx, xbp, xxa,
	    pot, xxf, xxp, xx0a, xx0f, xx0p, duma, dumb, rhoa, rhob, epsc,
	    rho16, zeta, rhot, coef1, coef2, coef3, xbq2a, xbq2f, zeta2,
	    zeta3, zeta4, dpot1, dpot2, xbq2p, epsca, depsc, epscf, epscp,
	    gzeta, coef1a, coef2a, coef3a, coef1f, coef2f, coef3f, coef1p,
	    coef2p, coef3p, bracka, brackf, depsca, brackp, depscf, depscp,
	    dgzeta, facvwn, curlya, curlyf, curlyp;




	/*                FACR = (THREE/(FOUR*PI))**SIXTH */
	/*                FACH = FOUR/((TWO**THIRD-ONE)*ANINE) */

	/*        THE A,B,C,X0 PARAMETERS OF THE FIT ARE AS FOLLOWS. */
	/*        P,F,A MEAN PARAMAGNETIC, FERROMAGNETIC, AND STIFFNESS TERMS. */
	/*        NOTE UNITS OF THE -A- CONSTANTS ARE RYDBERGS, SO DIVIDE BY 2. */



	/*     VWN CORRELATION FUNCTIONAL - FORMULA V */
	/*     S.H.VOSKO, L.WILK, M.NUSAIR  CAN.J.PHYS. 58, 1200-11(1980) */
	/*     THE NECESSARY FORMULAE ARE CLEARLY SPELLED OUT IN THE APPENDIX OF */
	/*     B.G.JOHNSON, P.M.W.GILL, J.A.POPLE  J.CHEM.PHYS. 98, 5612-26(1993) */
	/*     THE PARAMETERS APPEARING ABOVE CAN BE FOUND IN THE MOLPRO MANUAL, */
	/*     ALTHOUGH A FEW ARE VISIBLE IN THE VWN PAPER ITSELF. */

	/*     ON ENTRY, */
	/*     ROA,B   .... SPIN DENSITIES */
	/*     FTOTWT  .... VOLUME OF THIS GRID POINT */

	/*     ON EXIT, */
	/*     EVWN    .... CORRELATION ENERGY OF VWN5 */
	/*     VXCA,B  .... FIRST DERIVATIVES OF THE CORRELATION FUNCTIONAL */

	if (rhfgvb) {
	    facvwn = scale;
	    rhoa = roa;
	    rhob = rob;
	    rhot = rhoa + rhob;
	    rho16 = pow(rhot, c_b2);

	    q = sqrt(37.846991046399999);
	    xx0 = 12.5549141492;
	    coef1 = .031090699999999999;
	    coef2 = .012165997700463841 / xx0;
	    coef3 =
		7.4548800000000002 / q*.031090699999999999 *
		(12.924179199599999 / xx0);

	    x = .78762331789974325053 / rho16;
	    xx = x*x + x*3.72744 + 12.9352;
	    t1 = log(x*x / xx);
	    t2 = log((x + .10498)*(x + .10498) / xx);
	    t3 = atan(q / (x*2. + 3.72744));
	    epsc = coef1*t1 + coef2*t2 + coef3*t3;
	    evwn += facvwn*epsc*rhot*ftotwt;

	    depsc = (1. / x - x / xx*(3.72744 / (x + .10498) + 1.)) *
		.031090699999999999;
	    duma = epsc - x*depsc / 3.;
	    dumb = 0.;
	    vxca += facvwn*duma;
	    vxcb += dumb;
	}

	if (urohf) {
	    facvwn = scale;
	    rhoa = roa;
	    rhob = rob;
	    rhot = rhoa + rhob;
	    rho16 = pow(rhot, c_b2);

	    qp = sqrt(37.846991046399999);
	    qf = sqrt(22.381669423600009);
	    qa = sqrt(50.738680655099998);
	    xx0p = 12.5549141492;
	    xx0f = 15.868788500000001;
	    xx0a = 12.99914055888256;

	    /*          THE THREE POTENTIALS */

	    coef1p = .031090699999999999;
	    coef1f = .015545349999999999;
	    coef1a = -.01688686;
	    coef2p = .012165997700463841 / xx0p;
	    coef2f = .035670927515274994 / xx0f;
	    coef2a = -9.0886490370167686e-5 / xx0a;
	    coef3p =
		7.4548800000000002 / qp*.031090699999999999 *
		(12.924179199599999 / xx0p);
	    coef3f =
		14.120839999999999 / qf*.015545349999999999 *
		(17.952175 / xx0f);
	    coef3a =
		2.26214 / qa*-.01688686*(13.00447735762944 / xx0a);

	    x = .78762331789974325053 / rho16;
	    xxp = x*x + x*3.72744 + 12.9352;
	    xxf = x*x + x*7.06042 + 18.0578;
	    xxa = x*x + x*1.13107 + 13.0045;
	    t1p = log(x*x / xxp);
	    t1f = log(x*x / xxf);
	    t1a = log(x*x / xxa);
	    t2p = log((x + .10498)*(x + .10498) / xxp);
	    t2f = log((x + .325)*(x + .325) / xxf);
	    t2a = log((x + .0047584)*(x + .0047584) / xxa);
	    t3p = atan(qp / (x*2. + 3.72744));
	    t3f = atan(qf / (x*2. + 7.06042));
	    t3a = atan(qa / (x*2. + 1.13107));

	    epscp = coef1p*t1p + coef2p*t2p + coef3p*t3p;
	    epscf = coef1f*t1f + coef2f*t2f + coef3f*t3f;
	    epsca = coef1a*t1a + coef2a*t2a + coef3a*t3a;

	    /*          DERIVATIVE OF THE THREE POTENTIALS */

	    xbp = x*2. + 3.72744;
	    xbf = x*2. + 7.06042;
	    xba = x*2. + 1.13107;
	    xbq2p = xbp*xbp + qp*qp;
	    xbq2f = xbf*xbf + qf*qf;
	    xbq2a = xba*xba + qa*qa;
	    brackp = 2. / (x + .10498) - xbp / xxp - 14.06992 / xbq2p;
	    brackf =
		2. / (x + .325) - xbf / xxf - 25.641679999999997 / xbq2f;
	    bracka =
		2. / (x + .0047584) - xba / xxa -
		4.4862127999999997 / xbq2a;
	    curlyp = 2. / x - xbp / xxp - 14.90976 / xbq2p;
	    curlyf = 2. / x - xbf / xxf - 28.241679999999999 / xbq2f;
	    curlya = 2. / x - xba / xxa - 4.5242800000000001 / xbq2a;
	    depscp = (curlyp - brackp*-.39130665120000002 / xx0p) *
		.031090699999999999;
	    depscf = (curlyf - brackf*-2.2946365000000002 / xx0f) *
		.015545349999999999;
	    depsca = (curlya - bracka*-.0053820834880000008 / xx0a) *
		-.01688686;

	    /*          SOME AUXILIARY COMBINATIONS */

	    hx = (epscf - epscp)*1.70992093416136561756 / epsca - 1.;
	    dhx = (depscf - depscp - (epscf - epscp)*depsca / epsca) *
		1.70992093416136561756 / epsca;

	    zeta = (rhoa - rhob) / rhot;
	    zeta2 = zeta*zeta;
	    zeta3 = zeta2*zeta;
	    zeta4 = zeta2*zeta2;
	    d__1 = zeta + 1.;
	    d__2 = 1. - zeta;
	    gzeta = (pow(d__1, c_b4) + pow(d__2, c_b4) - 2.)*1.125;
	    d__1 = zeta + 1.;
	    d__2 = 1. - zeta;
	    dgzeta = (pow(d__1, c_b6) - pow(d__2, c_b6))*1.5;

	    /*          TOTAL POTENTIAL, AND CONTRIBUTION TO ENERGY */

	    pot = epscp + epsca*gzeta*(hx*zeta4 + 1.);
	    evwn += facvwn*pot*rhot*ftotwt;

	    /*          CONTRIBUTION TO FUNCTIONAL DERIVATIVE */

	    dpot1 =
		-(x / 6.)*(depscp + depsca*gzeta*(hx*zeta4 + 1.) +
			     epsca*gzeta*dhx*zeta4);
	    dpot2 =
		epsca*(dgzeta*(hx*zeta4 + 1.) +
			 gzeta*4.*hx*zeta3);
	    duma = pot + (dpot1 + dpot2*(1. - zeta));
	    dumb = pot + (dpot1 - dpot2*(zeta + 1.));

	    vxca += facvwn*duma;
	    vxcb += facvwn*dumb;
	}

	return 0;
    }
};
