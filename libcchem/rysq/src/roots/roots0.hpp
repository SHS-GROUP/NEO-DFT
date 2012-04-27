#include <cmath>

// Polynomial Coefficients - Parallelize root calculations

namespace roots0 {

    static const TYPE__ PIE4 = 0.785398163397448;

    // BOOST_GPU_ENABLED
    // inline TYPE__ Ghondo_Rysq_roots0(TYPE__ X) {
  
    // 	//#if defined (USE_GHONDO)
  
    // 	if (X > 33.) {
    // 	    return Ghondo_roots0_x33toInf(X);
    // 	}
    // 	if (X > 15.) {
    // 	    return Ghondo_roots0_x15to33(X);
    // 	}
    // 	if (X > 10.) {
    // 	    return Ghondo_roots0_x10to15(X);
    // 	}
    // 	if (X > 5.) {
    // 	    return Ghondo_roots0_x5to10(X);
    // 	}
    // 	if (X > 3.) {
    // 	    return Ghondo_roots0_x3to5(X);
    // 	}
    // 	if (X > 1.) {
    // 	    return Ghondo_roots0_x1to3(X);
    // 	}
    // 	if (X > 3e-7) {
    // 	    return Ghondo_roots0_x0to1(X);
    // 	}
    // 	return (1. - X / 3.);
    // }  

    //#else // Explicit Polynomial Expansion
    BOOST_GPU_ENABLED
    inline void evaluate(TYPE__ X, TYPE__ (&t2)[0], TYPE__ (&W_)[1]) {
    
	const TYPE__ C15_33[] = { -6.0156581186480998977e-5,
				  -4.9695241464489997574e-1,
				  1.9623264149430000303e-1 };


	const TYPE__ C10_15[] = { -2.1916512131607000959e-5,
				  -4.9893752514047001734e-1,
				  2.2991849164985000975e-1,
				  -1.8784686463511998666e-1 };

	const TYPE__ C5_10[] = { -3.1501078774085e-6,
				 -.49984072848436,
				 .24645596956002,
				 -.32883030418398,
				 .53689283271887,
				 -.69955602298985,
				 .46897511375022};

	const TYPE__ C3_5[] = { .0528406320615584,
				-.0175257821619926,
				.00433207949514611,
				-8.81947375894379e-4,
				1.526537461148e-4,
				-2.290098979647e-5,
				3.022556449731e-6,
				-3.553558319675e-7,
				3.760256799971e-8,
				-3.614965656163e-9,
				3.24031041623823e-10,
				-2.62453564772299e-11};

	const TYPE__ C1_3[] = {.115702180856167,
			       -.0529428148329736,
			       .0161723488664661,
			       -.00379490003707156,
			       7.24888732052332e-4,
			       -1.16740298039895e-4,
			       1.62429321438911e-5,
			       -1.98850171329371e-6,
			       2.17216556336318e-7,
			       -2.14234468198419e-8,
			       1.96215250865776e-9,
			       -1.61702782425558e-10};

	const TYPE__ C0_1[] = {.333333333333318,
			       -.199999999997023,
			       .0714285713298222,
			       -.0185185172458485,
			       .00378787044215009,
			       -6.40994113129432e-4,
			       9.25197374512647e-5,
			       -1.15662609053481e-5,
			       1.21222603512827e-6,
			       -8.36313918003957e-8};


	TYPE__ &W = W_[0];
	W = 0.0;

	if(X > 5.) {
	    TYPE__ x1 = 1.0/X;
	    TYPE__ x1i = 1.0;
	    if(X > 33.) { }
	    else if(X > 15.) {
		const int lsize = sizeof(C15_33)/sizeof(TYPE__);
		for(int i = 0; i < lsize; ++i) {
		    W += x1i*C15_33[i];
		    x1i *= x1;
		}
	    }
	    else if(X > 10.) {
		const int lsize = sizeof(C10_15)/sizeof(TYPE__);
		for(int i = 0; i < lsize; ++i) {
		    W += x1i*C10_15[i];
		    x1i *= x1;
		}
	    
	    }
	    else {
		const int lsize = sizeof(C5_10)/sizeof(TYPE__);
		for(int i = 0; i < lsize; ++i) {
		    W += x1i*C5_10[i];
		    x1i *= x1;
		}
	    }
	
	    W *= exp(-X);
	    W += sqrt(PIE4*x1);
	    return;
	}
	else if(X > 3e-7) {
	    TYPE__ yi = 1.0;
	    TYPE__ y = 0.0;
	    if(X > 3.0) {
		y =  X - 4.;
		int lsize = sizeof(C3_5)/sizeof(TYPE__);
		for(int i = 0; i < lsize;  ++i) {
		    W += C3_5[i]*yi;
		    yi *= y;
		}
	    }
	    else if(X > 1.0) {
		y = X - 2.;
		int lsize = sizeof(C1_3)/sizeof(TYPE__);
		for(int i = 0; i < lsize;  ++i) {
		    W += C1_3[i]*yi;
		    yi *= y;
		}
	    }
	    else {
		y = X;
		int lsize = sizeof(C0_1)/sizeof(TYPE__);
		for(int i = 0; i < lsize; ++i) {
		    W += C0_1[i]*yi;
		    yi *= y;
		}
	    }
	    W *= (X + X);
	    W += exp(-X);
	    return;
	}

	W = (1.0 - X*(1.0/3.0));

    }

} // namespace roots0

