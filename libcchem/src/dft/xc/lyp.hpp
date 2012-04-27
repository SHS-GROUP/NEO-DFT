#include "dft/xc/kernel.hpp"
#include "dft/xc/gamess.hpp"

struct lyp : needs_gradient<1> {

    const double scale_;
    explicit lyp(double scale) : scale_(scale) {} 

    void operator()(const Point::Density &density, Point::XC &xc) {
	gamess::correlation(&lypcf_, scale_, density, xc);
    }

    static 
    int lypcf_(const double scale,
	       const bool rhfgvb, const bool urohf, const int nb,
	       const double roa, const double rob, const double grdaa,
	       const double grdbb, const double grdab, const double gradxa,
	       const double gradya, const double gradza,
	       const double gradxb, const double gradyb,
	       const double gradzb, const double ftotwt, double &elyp,
	       double &duma, double &dumxa, double &dumya, double &dumza,
	       double &dumb, double &dumxb, double &dumyb, double &dumzb) {

	const double c_b2 = 2.;
	const double c_b3 = 3.6666666666666665;
	const double c_b4 = .66666666666666663;
	const double c_b5 = .33333333333333331;

	double d__1, d__2;

	/* Local variables */
	double u, w, t1, t2, ra, rb, ua, ub, xc, xd, rx, t2a, t2b, faa,
	    gaa, gab, gbb, fbb, fab, dda, ddb, raa, rbb, rab, waa, wab,
	    wbb, ddaa, ddab, ddbb, face, facp, rhoa, rhob, rho13, rho43,
	    zeta, rhot, dt1da, dt2da, dt1db, dt2db, omega, delta, rho13a,
	    rho13b, rho83a, alypa, alypb, alypc, alypd, rho83b, elypi,
	    d2daaa, d2dabb, d2daab, d2dbaa, d2dbbb, d2dbab, dfaada, dfabda,
	    dfbbda, dfaadb, dfbbdb, dfabdb, deltap, omegap;

	/* ************************************ */
	/*     LYP CORRELATION FUNCTIONAL */

	/*     C. LEE, W. YANG AND R. G. PARR */
	/*     PHYSICAL REVIEW B 37, 785-789 (1988) */

	/*     ----- */
	/*     INPUT */
	/*     ----- */

	/*     ROA,B   .... SPIN DENSITIES */
	/*     GRAD*  .... DENSITY GRADIENTS */
	/*     FTOTWT  .... VOLUME OF THIS GRID */

	/*     ------ */
	/*     OUTPUT */
	/*     ------ */

	/*     ELYP    ....  CORRELATION ENERGY OF LYP */
	/*     DUM**  ....  FIRST DERIVATIVES OF THE EXCHANGE FUNCTINALS */
	/* ************************************ */

	/* Computing 2nd power */
	d__2 = PI;
	d__1 = d__2 * d__2 * 3.;
	face = pow(c_b2, c_b3) * .29999999999999999 * pow(d__1, c_b4);
	facp = scale * 1.;
	alypa = .04918;
	alypb = .132;
	alypc = .2533;
	alypd = .349;
	/* ************************************ */
	if (rhfgvb) {
	    /* ************************************ */
	    rhoa = roa;
	    rhob = rob;
	    gaa = grdaa;
	    gbb = grdbb;
	    gab = grdab;
	    rhot = rhoa + rhob;
	    rho13 = pow(rhot, c_b5);
	    rho13a = pow(rhoa, c_b5);
	    rho13b = pow(rhob, c_b5);
	    /* Computing 4th power */
	    d__1 = rho13, d__1 *= d__1;
	    rho43 = d__1 * d__1;
	    /* Computing 8th power */
	    d__1 = rho13a, d__1 *= d__1, d__1 *= d__1;
	    rho83a = d__1 * d__1;
	    /* Computing 8th power */
	    d__1 = rho13b, d__1 *= d__1, d__1 *= d__1;
	    rho83b = d__1 * d__1;
	    xc = alypc / rho13;
	    xd = alypd / rho13;
	    rx = xd / (xd + 1.);
	    delta = xc + rx;
	    omega = alypa * alypb / alypd * rx * exp(-xc) / rho43;
	    deltap = (-delta + rx * rx) / 3. / rhot;
	    omegap = (delta - 5.) / 3. / rhot;
	    ra = rhoa / rhot;
	    rb = rhob / rhot;
	    raa = ra * ra;
	    rbb = rb * rb;
	    rab = ra * rb;
	    /* ************************************ */
	    faa = rab / 9 * (1. - delta * 3. - (delta - 11) * ra);
	    fbb = rab / 9 * (1. - delta * 3. - (delta - 11) * rb);
	    fab = rab / 9 * (47 - delta * 7);
	    ddaa = -omega * (faa - rbb);
	    ddbb = -omega * (fbb - raa);
	    ddab = -omega * (fab - 1.3333333333333333);
	    t1 = alypa * -4. / alypd * rx * rab * rho43;
	    t2a = -face * omega * rab * rho83a;
	    t2b = -face * omega * rab * rho83b;
	    t2 = t2a + t2b;
	    elypi = t1 + t2 + ddaa * gaa + ddbb * gbb + ddab * gab;
	    /* ************************************ */
	    dfaada =
		-(rab / 9) * ((ra + 3.) * deltap +
			      (delta - 11) * rb / rhot);
	    dfbbda =
		-(rab / 9) * ((rb + 3.) * deltap -
			      (delta - 11) * rb / rhot);
	    dfabda = -(rab / 9) * 7 * deltap;
	    d2daaa = omegap * ddaa - omega * (dfaada + rbb * (2. / rhot));
	    d2dabb = omegap * ddbb - omega * (dfbbda - rab * (2. / rhot));
	    d2daab = omegap * ddab - omega * dfabda;
	    dt1da = t1 / rhot * (rx / 3. + 1.);
	    dt2da = omegap * (t2a + t2b) + t2a * 2.6666666666666665 / rhoa;
	    dda =
		dt1da + dt2da + d2daaa * gaa + d2dabb * gbb + d2daab * gab;
	    /* ************************************ */
	    elyp += facp * elypi * ftotwt;
	    u = dda;
	    u *= facp;
	    duma += u;
	    w = ddaa + ddab * .5;
	    w *= facp;
	    dumxa += w * 2. * gradxa;
	    dumya += w * 2. * gradya;
	    dumza += w * 2. * gradza;
	    /* ************************************ */
	} else if (urohf) {
	    /* ************************************ */
	    gaa = grdaa;
	    gbb = grdbb;
	    gab = grdab;
	    rhoa = roa;
	    rhob = rob;
	    rhot = rhoa + rhob;
	    rho13 = pow(rhot, c_b5);
	    rho13a = pow(rhoa, c_b5);
	    rho13b = pow(rhob, c_b5);
	    /* Computing 4th power */
	    d__1 = rho13, d__1 *= d__1;
	    rho43 = d__1 * d__1;
	    /* Computing 8th power */
	    d__1 = rho13a, d__1 *= d__1, d__1 *= d__1;
	    rho83a = d__1 * d__1;
	    /* Computing 8th power */
	    d__1 = rho13b, d__1 *= d__1, d__1 *= d__1;
	    rho83b = d__1 * d__1;
	    xc = alypc / rho13;
	    xd = alypd / rho13;
	    rx = xd / (xd + 1.);
	    delta = xc + rx;
	    omega = alypa * alypb / alypd * rx * exp(-xc) / rho43;
	    deltap = (-delta + rx * rx) / 3. / rhot;
	    omegap = (delta - 5.) / 3. / rhot;
	    ra = rhoa / rhot;
	    rb = rhob / rhot;
	    raa = ra * ra;
	    rbb = rb * rb;
	    rab = ra * rb;
	    zeta = ra - rb;
	    /* ************************************ */
	    faa = rab / 9 * (1. - delta * 3. - (delta - 11) * ra);
	    fbb = rab / 9 * (1. - delta * 3. - (delta - 11) * rb);
	    fab = rab / 9 * (47 - delta * 7);
	    ddaa = -omega * (faa - rbb);
	    ddbb = -omega * (fbb - raa);
	    ddab = -omega * (fab - 1.3333333333333333);
	    t1 = alypa * -4. / alypd * rx * rab * rho43;
	    t2a = -face * omega * rab * rho83a;
	    t2b = -face * omega * rab * rho83b;
	    t2 = t2a + t2b;
	    elypi = t1 + t2 + ddaa * gaa + ddbb * gbb + ddab * gab;
	    /* ************************************ */
	    dfaada =
		-zeta * faa / rhoa - rab / 9 * ((ra + 3.) * deltap +
						(delta - 11) * rb / rhot);
	    dfbbda =
		-zeta * fbb / rhoa - rab / 9 * ((rb + 3.) * deltap -
						(delta - 11) * rb / rhot);
	    dfabda = -zeta * fab / rhoa - rab / 9 * 7 * deltap;
	    if (nb > 0) {
		dfaadb =
		    zeta * faa / rhob - rab / 9 * ((ra + 3.) * deltap -
						   (delta -
						    11) * ra / rhot);
		dfbbdb =
		    zeta * fbb / rhob - rab / 9 * ((rb + 3.) * deltap +
						   (delta -
						    11) * ra / rhot);
		dfabdb = zeta * fab / rhob - rab / 9 * 7 * deltap;
	    } else {
		dfaadb = 0.;
		dfbbdb = 0.;
		dfabdb = 0.;
	    }
	    /* ************************************ */
	    d2daaa = omegap * ddaa - omega * (dfaada + rbb * (2. / rhot));
	    d2dbaa = omegap * ddaa - omega * (dfaadb - rab * (2. / rhot));
	    d2dabb = omegap * ddbb - omega * (dfbbda - rab * (2. / rhot));
	    d2dbbb = omegap * ddbb - omega * (dfbbdb + raa * (2. / rhot));
	    d2daab = omegap * ddab - omega * dfabda;
	    d2dbab = omegap * ddab - omega * dfabdb;
	    /* ************************************ */
	    dt1da = t1 / rhot * (rx / 3. + rb / ra);
	    dt2da =
		(omegap - zeta / rhoa) * t2b + (omegap +
						(2.6666666666666665 -
						 zeta) / rhoa) * t2a;
	    if (nb > 0) {
		dt1db = t1 / rhot * (rx / 3. + ra / rb);
		dt2db = (omegap + zeta / rhob) * t2a + (omegap + (zeta +
								  2.6666666666666665)
							/ rhob) * t2b;
	    } else {
		dt1db = 0.;
		dt2db = 0.;
	    }
	    dda =
		dt1da + dt2da + d2daaa * gaa + d2dabb * gbb + d2daab * gab;
	    ddb =
		dt1db + dt2db + d2dbaa * gaa + d2dbbb * gbb + d2dbab * gab;

	    elyp += facp * elypi * ftotwt;
	    ua = dda;
	    ua *= facp;
	    duma += ua;
	    ub = ddb;
	    ub *= facp;
	    dumb += ub;
	    waa = ddaa * 2.;
	    waa *= facp;
	    wbb = ddbb * 2.;
	    wbb *= facp;
	    wab = ddab;
	    wab *= facp;
	    dumxa = dumxa + waa * gradxa + wab * gradxb;
	    dumya = dumya + waa * gradya + wab * gradyb;
	    dumza = dumza + waa * gradza + wab * gradzb;
	    dumxb = dumxb + wbb * gradxb + wab * gradxa;
	    dumyb = dumyb + wbb * gradyb + wab * gradya;
	    dumzb = dumzb + wbb * gradzb + wab * gradza;
	}
	return 0;
    }

};
