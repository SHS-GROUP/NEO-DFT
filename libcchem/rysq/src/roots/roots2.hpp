#include "roots/evaluate.hpp"

namespace roots2 {

    static const TYPE__ PIE4 = 0.785398163397448;

    static const int N = 2;

    constant__
    TYPE__ R2[] = { 0.275255128608411,  2.72474487139158 };

    constant__
    TYPE__ W2[] = { 1.0,  0.0917517095361369 };

    namespace x40toInf {
	template<class A, typename T>
	BOOST_GPU_ENABLED
	void evaluate(TYPE__ X, A &R, A &W, TYPE__ W0,
		      const T &thread) {
	    evaluate_inf( R, X, R2, W, W0,  W2, thread);
	}
    }

    namespace x33to40 {
	constant__
	TYPE__ CR[2][2] = {{ -0.87894730749888, 10.9243702330261 },
			   { -9.28903924275977, 81.0642367843811 }};
	constant__
	TYPE__ CW[2][2] = {{ 0.0, 0.0 },
			   { 14.468573893084, -77.9250653461045 }};
	template<class A, typename T>
	BOOST_GPU_ENABLED 
	void evaluate(TYPE__ X, A &R, A &W,
		      TYPE__ W0, TYPE__ e, const T &thread) {
	    polynomial::evaluate2(X, R, W, CR, CW, thread);
	    scale_add(e, R, W, X, R2, W0, W2, thread);
	    if (master(thread)) W[0] -= W[1];//sum<N-1>(W+1);
	}
    }

    namespace x15to33 {

	constant__
	TYPE__ CR1[2][3] = {{ -0.0231080873898939, 9.21005186542857, -47.5742064274859 },
			    { 2.98011277766958, -86.6891724287962, -135.831002139173 }};

	constant__
	TYPE__ CR2[2][4] = {{ -0.137292644149838, -0.0171984023644904, 1.76003409708332E-4, -1.14906395546354E-6 },
			    { -4.02886174850252, -0.0971850973831558, 3.64921633404158E-4, 0.0 }};

	constant__
	TYPE__ CF[3] = { -1.0 - 6.0156581186481E-5, -0.496952414644900, 0.196232641494300 };

	template<class A, typename T>
	BOOST_GPU_ENABLED
	void evaluate(TYPE__ X, A &R, TYPE__ &f,
		      TYPE__ W0, TYPE__ e, const T &thread) {
	    polynomial::inverse(R, CR1, X, thread);
	    add_evaluate_polynomial1(R, CR2, X, thread);
	    scale<N>(e, R, thread);
	    add_evaluate_inf(R, R2, X, thread);
	    if (master(thread))
		f = (polynomial::inverse(CF, X)*e + W0)/(2*X);
	}
    }

    namespace x10to15 {
	constant__
	TYPE__ CF[4] = { -1.0 - 2.1916512131607E-5, -0.498937525140470, 0.229918491649850, -0.187846864635120 };
	constant__
	TYPE__ CR2[2][4] = {{ 1.25705571069895, -0.0673760231824074, 0.00119483054115173, -1.01041157064226E-5 },
			    { -4.22216483306320, -0.0934976436343509, 3.39024225137123E-4, 0.0 }};

	constant__
	TYPE__ CR1[2][5] = {{ -23.8570496490846, 264.536689959503, -1708.07677109425, 5910.05939591842, -8576.09422987199 },
			    { 8.00839033297501, -156.184800325063, 339.891508992661, -1049.99071905664, -2084.57050986847 }};


	template<class A, typename T>
	BOOST_GPU_ENABLED
	void evaluate(TYPE__ X, A &R, TYPE__ &f, TYPE__ W0, TYPE__ e,
		      const T &thread) {
      	    polynomial::inverse(R, CR1, X, thread);
	    add_evaluate_polynomial1(R, CR2, X, thread);
	    scale<N>(e, R, thread);
	    add_evaluate_inf(R, R2, X, thread);
	    if (master(thread))
		f = (polynomial::inverse(CF, X)*e + W0)/(2*X);
	}
    }

    namespace x5to10 {
	constant__
	TYPE__ CR[2][15] = {{ 0.0372458587837249, -0.00460448960876139, 4.74791204651451E-4, -3.76953869614706E-5, 2.283485114279E-6, -1.497500187053E-7, 1.71382882399E-8, -1.775564159743E-9, 6.628721099436E-11, 7.802544789997E-12, -7.719300212748E-13, -7.064522786879E-14, 1.3583196188E-14, 2.38198922570405E-16, -1.43632730148572E-16 },
			    { 0.544765245686790, -0.0953478510453887, 0.0135094294917224, -0.00144614664924989, 1.0323664788832E-4, -2.4995220023291E-6, -4.383376014528E-7, 6.128153450169E-8, -2.624183464275E-9, -2.222722579924E-10, 4.190559455515E-11, -2.224334349799E-12, -1.36113510175724E-13, 2.487916227989E-14, 0.0 }};
	constant__
	TYPE__ CF[7] = { -1.0 - 3.1501078774085E-6, -0.499840728484360, 0.246455969560020, -0.328830304183980, 0.536892832718870, -0.699556022989850, 0.468975113750220 };

	template<class A, typename T>
	BOOST_GPU_ENABLED
	void evaluate(TYPE__ X, A &R, TYPE__ &f, TYPE__ W0, TYPE__ e,
		      const T &thread) {
	    polynomial::evaluate(R, TYPE__(X - 7.5), CR, thread);
	    if (master(thread))
		f = (polynomial::inverse(CF, X)*e + W0)/(2*X);
	}
    }


    namespace x3to5 {
	constant__
	TYPE__ CF[12] = { 0.0528406320615584, -0.0175257821619926, 0.00433207949514611, -8.81947375894379E-4, 1.526537461148E-4, -2.290098979647E-5, 3.022556449731E-6, -3.553558319675E-7, 3.760256799971E-8, -3.614965656163E-9, 3.24031041623823E-10, -2.62453564772299E-11 };
	constant__
	TYPE__ CR[2][11] = {{ 0.0612589701086408, -0.00989572595720351, 0.00116210907653515, -1.12608004981982E-4, 1.08484384385679E-5, -9.76085576741771E-7, 5.93066856324744E-8, -1.73508862390291E-9, 7.10910223886747E-11, -4.11560117487296E-12, 0.0 },
			    { 1.12155283108289, -0.260417417692375, 0.0356550690684281, -0.00249327728643089, -2.04685420150802E-5, 1.85882653064034E-5, -7.017002532106E-7, -1.497986283037E-7, 1.60349804524E-8, 5.44072475994123E-10, -1.80555625241001E-10 }};
	template<class A, typename T>
	BOOST_GPU_ENABLED
	void evaluate(TYPE__ X, A &R, TYPE__ &f, TYPE__ W0, TYPE__ e,
		      const T &thread) {
	    polynomial::evaluate(R, TYPE__(X - 4.0), CR, thread);
	    if (master(thread))
		f = evaluate_polynomial(CF, (X - 4));
	}
    }

    namespace x1to3 {
	constant__
	TYPE__ CF[12] = { 0.115702180856167, -0.0529428148329736, 0.0161723488664661, -0.00379490003707156, 7.24888732052332E-4, -1.16740298039895E-4, 1.62429321438911E-5, -1.98850171329371E-6, 2.17216556336318E-7, -2.14234468198419E-8, 1.96215250865776E-9, -1.61702782425558E-10 };
	constant__
	TYPE__ CR[2][11] = {{ 0.0868085688285261, -0.0163329339286794, 0.00219173220020161, -2.49018321709815E-4, 2.47191693238413E-5, -1.85306035634293E-6, 8.47225338838E-8, -3.846389873308E-10, -5.152207846962E-10, 8.4741706477627E-11, -6.36859636616415E-12 },
			    { 1.80400974537950, -0.430761584997596, 0.0485928174507904, -0.00158064140671893, -1.93160765581969E-4, 9.76783813082564E-6, 2.247389642339E-6, -1.725838516261E-7, -1.878920917404E-8, 2.07111465297976E-9, 1.45331350488343E-10 }};
	template<class A, typename T>
	BOOST_GPU_ENABLED
	void evaluate(TYPE__ X, A &R, TYPE__ &f, TYPE__ W0, TYPE__ e,
		      const T &thread) {
	    polynomial::evaluate(R, TYPE__(X - 2.0), CR, thread);
	    if (master(thread))
		f = evaluate_polynomial(CF, (X - 2));
	}
    }

    namespace x0to1 {
	constant__
	TYPE__ CF[10] = { 0.333333333333318, -0.199999999997023, 0.0714285713298222, -0.0185185172458485, 0.00378787044215009, -6.40994113129432E-4, 9.25197374512647E-5, -1.15662609053481E-5, 1.21222603512827E-6, -8.36313918003957E-8 };
	constant__
	TYPE__ CR[2][9] = {{ 0.130693606237085, -0.0290430236084697, 0.00444654947116579, -5.33184749432408E-4, 4.743292959463E-5, -2.447252174587E-6, -4.558315364581E-8, 2.49173650389842E-8, -2.35234358048491E-9 },
			   { 2.86930639376289, -0.637623643056745, 0.0532735082098139, -5.88154362858038E-5, -1.345693393936E-4, -2.066168802076E-5, 1.83536773631E-6, 2.36809910635906E-7, -2.4740490232917E-8 }};
	template<class A, typename T>
	BOOST_GPU_ENABLED
	void evaluate(TYPE__ X, A &R, TYPE__ &f, TYPE__ W0, TYPE__ e,
		      const T &thread) {
	    polynomial::evaluate(R, X, CR, thread);
	    if (master(thread))
		f = evaluate_polynomial(CF, X);
	}
    }

    namespace x0 {
	constant__
	TYPE__ CR[2][2] = {{ 0.130693606237085, -0.0290430236082028 },
			   { 2.86930639376291, -0.637623643058102 }};
	constant__
	TYPE__ CW[2][2] =  {{ 0.652145154862545, -0.122713621927067 },
			    { 0.347854845137453, -0.210619711404725 }};
	template<class A, typename T>
	BOOST_GPU_ENABLED
	void evaluate(TYPE__ X, A &R, A &W,
		      const T &thread) {
	    polynomial::evaluate(R, X, CR, thread);
	    polynomial::evaluate(W, X, CW, thread);
	}
    }

    template<class A, typename T>
    BOOST_GPU_ENABLED
    void evaluate_(TYPE__ X, A &R, A &W, const T &thread) {
    
	const TYPE__ W0 = sqrt(PIE4/X); 
	TYPE__ e = exp(-X);

	if (X > 40.0) {
	    x40toInf::evaluate(X, R, W, W0, thread);
	    goto end;
	}
	else if (X > 33.0) {
	    x33to40::evaluate(X, R, W, W0, e, thread);
	    goto end;
	}
	else if (X > 15.0) {
	    x15to33::evaluate(X, R, W[0], W0, e, thread);
	}  
	else if (X > 10.0) {
	    x10to15::evaluate(X, R, W[0], W0, e, thread);
	}
	else if (X > 5.0) {
	    x5to10::evaluate(X, R, W[0], W0, e, thread);
	}
	else if (X > 3.0) {
	    x3to5::evaluate(X, R, W[0], W0, e, thread);
	}
	else if (X > 1.0) {
	    x1to3::evaluate(X, R, W[0], W0, e, thread);
	}
	else if (X > 3.0e-7) {
	    x0to1::evaluate(X, R, W[0], W0, e, thread);
	}
	else {
	    x0::evaluate(X, R, W, thread);
	    goto end;
	}

	if (master(thread)) {
	    TYPE__ f = W[0];
	    W[0] = 2*X*f + e;
	    W[1] = ((f - W[0])*R[0] + f)*(R[1] + 1.0)/(R[1] - R[0]);
	    W[0] -= W[1];
	}


    end:

	change_variable<N>(R, thread);
	// for(int i = 0; i < 2; ++i) {
	// 	R[i] = R[i]/(1.0 + R[i]);
	// }

	return;
    }


    BOOST_GPU_ENABLED
    inline void evaluate(TYPE__ X, TYPE__ (&R)[2], TYPE__ (&W)[2]) {
	evaluate_(X, R, W, serial_tag());
    }

    template<typename T>
    BOOST_GPU_ENABLED
    inline void evaluate(TYPE__ X, TYPE__ *R, TYPE__ *W, const T &thread) {
	evaluate_(X, R, W, thread);
    }

} // namespace roots2
