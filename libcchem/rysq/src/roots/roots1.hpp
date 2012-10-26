#include <math.h>

namespace roots1 {

    static const TYPE__ PIE4 = 0.785398163397448;

    static const TYPE__ r12 = .275255128608411;
    static const TYPE__ pie4 = .785398163397448;
    static const TYPE__ r22 = 2.72474487139158;
    static const TYPE__ w22 = .0917517095361369;
    static const TYPE__ r13 = .190163509193487;
    static const TYPE__ r23 = 1.78449274854325;
    static const TYPE__ w23 = .177231492083829;
    static const TYPE__ r33 =  5.52534374226326;
    static const TYPE__ w33 = .00511156880411248;


    BOOST_GPU_ENABLED
    inline void x33toINF(TYPE__ X, TYPE__ *t2, TYPE__ *W) {
	W[0] = sqrt(pie4 / X);
	t2[0] = .5 / (X - .5);
    }

    BOOST_GPU_ENABLED
    inline void x15to33(TYPE__ X, TYPE__ *t2, TYPE__ *W) {
	TYPE__ e = exp(-X);
	W[0] = ((.1962326414943 / X - .4969524146449) / X - 
		6.0156581186481e-5) * e + sqrt(pie4 / X);
	TYPE__ f1 = (W[0] - e) / (X + X);
	t2[0] = f1 / (W[0] - f1);
    }

    BOOST_GPU_ENABLED
    inline void x10to15(TYPE__ X, TYPE__ *t2, TYPE__ *W) {
	TYPE__ e = exp(-X);
	W[0] = (((-.18784686463512 / X + .22991849164985) / X - 
		 .49893752514047) / X - 2.1916512131607e-5) * e + sqrt(pie4/X);
	TYPE__ f1 = (W[0] - e) / (X + X);
	t2[0] = f1 / (W[0] - f1);
    }

    BOOST_GPU_ENABLED
    inline void x5to10(TYPE__ X, TYPE__ *t2, TYPE__ *W) {
	TYPE__ e = exp(-X);
	W[0] = ((((((.46897511375022 / X - .69955602298985) / X + 
		    .53689283271887) / X - .32883030418398) / X + 
		  .24645596956002) / X - .49984072848436) / X - 
		3.1501078774085e-6) * e + sqrt(pie4 / X);
	TYPE__ f1 = (W[0] - e) / (X + X);
	t2[0] = f1 / (W[0] - f1);
    }

    BOOST_GPU_ENABLED
    inline void x3to5(TYPE__ X, TYPE__ *t2, TYPE__ *W) {
	TYPE__ y = X - 4.;
	TYPE__ f1 = ((((((((((y * -2.62453564772299e-11 + 3.24031041623823e-10) * y - 
			     3.614965656163e-9) * y + 3.760256799971e-8) * y - 
			   3.553558319675e-7) * y + 3.022556449731e-6) * y - 
			 2.290098979647e-5) * y + 1.526537461148e-4) * y - 
		       8.81947375894379e-4) * y + .00433207949514611) * y - 
		     .0175257821619926) * y + .0528406320615584;
	W[0] = (X + X) * f1 + exp(-X);
	t2[0] = f1 / (W[0] - f1);
    }

    BOOST_GPU_ENABLED
    inline void x1to3(TYPE__ X, TYPE__ *t2, TYPE__ *W) {
	TYPE__ y = X - 2.;
	TYPE__ f1 = ((((((((((y * -1.61702782425558e-10 + 1.96215250865776e-9) * y - 
			     2.14234468198419e-8) * y + 2.17216556336318e-7) * y - 
			   1.98850171329371e-6) * y + 1.62429321438911e-5) * y - 
			 1.16740298039895e-4) * y + 7.24888732052332e-4) * y - 
		       .00379490003707156) * y + .0161723488664661) * y - 
		     .0529428148329736) * y + .115702180856167;
	W[0] = (X + X) * f1 + exp(-X);
	t2[0] = f1 / (W[0] - f1);
    }

    BOOST_GPU_ENABLED
    inline void x0to1(TYPE__ X, TYPE__ *t2, TYPE__ *W) {
	TYPE__ f1 = ((((((((X * -8.36313918003957e-8 + 1.21222603512827e-6) * 
			   X - 1.15662609053481e-5) * X + 9.25197374512647e-5) 
			 * X - 6.40994113129432e-4) * X + .00378787044215009)
		       * X - .0185185172458485) * X + .0714285713298222) *
		     X - .199999999997023) * X + .333333333333318;
	W[0] = (X + X) * f1 + exp(-X);
	t2[0] = f1 / (W[0] - f1);
    }

    template<class A, typename T>
    BOOST_GPU_ENABLED
    void evaluate_(TYPE__ X, A &t2, A &W, const T &thread) {

	if (!master(thread)) return;

	if (X > 33.) {
	    x33toINF(X, t2, W);
	}
	else if(X > 15.) {
	    x15to33(X, t2, W);
	}
	else if (X > 10.) {
	    x10to15(X, t2, W);
	}
	else if (X > 5.) {
	    x5to10(X, t2, W);
	}
	else if (X > 3.) {
	    x3to5(X, t2, W);
	}
	else if (X > 1.) {
	    x1to3(X, t2, W);
	}
	else if (X > 3e-7) {
	    x0to1(X, t2, W);
	}
	/*     X IS APPROXIMATELY ZERO.         NROOTS=1 */
	else {
	    t2[0] = .5 - X / 5.;
	    W[0] = 1. - X / 3.;
	}

	for(int i = 0; i < 1; ++i) {
	    t2[i] = t2[i]/(1.0 + t2[i]);
	}
	return;
    }


    BOOST_GPU_ENABLED
    inline void evaluate(TYPE__ X, TYPE__ (&R)[1], TYPE__ (&W)[1]) {
	evaluate_(X, R, W, serial_tag());
    }

    template<typename T>
    BOOST_GPU_ENABLED
    inline void evaluate(TYPE__ X, TYPE__ *R, TYPE__ *W, const T &thread) {
	evaluate_(X, R, W, thread);
    }

} // namespace roots1

