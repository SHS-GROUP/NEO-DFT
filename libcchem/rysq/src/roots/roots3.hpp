#include "roots/evaluate.hpp"
#include <cstdlib>
#include <cmath>

namespace roots3 {

    static const TYPE__ PIE4 = 0.785398163397448;

    static const size_t N = 3;
    
    constant__
    TYPE__ R3[3] = { 0.190163509193487, 1.78449274854325, 5.52534374226326 };
    constant__
    TYPE__ W3[3] = { 1.0, 0.177231492083829, 0.00511156880411248 };

    namespace x0 {
	constant__
	TYPE__ R_CP[3][2] = {{ 0.0603769246832797, -0.00928875764357368 },
			     { 0.776823355931043, -0.119511285527878 },
			     { 6.66279971938567, -1.02504611068957 }};
	constant__
	TYPE__ W_CP[3][2] = {{ 0.467913934572691, -0.0564876917232519 },
			     { 0.360761573048137, -0.149077186455208 },
			     { 0.171324492379169, -0.127768455150979 }};

	BOOST_GPU_ENABLED
	static void evaluate(TYPE__ X, TYPE__ (&R)[3], TYPE__ (&W)[3],
			     const serial_tag& = serial_tag()) {
	    polynomial::evaluate2(X, R, W, R_CP , W_CP );
	}

	template<class A, typename  T>
	BOOST_GPU_ENABLED
	static void evaluate(TYPE__ X, A &R, A &W,
			     const T &thread) {
	    polynomial::evaluate2(X, R, W, R_CP , W_CP, thread);
	}

    }

    namespace x0to1 {
	constant__
	TYPE__  R_CP[3][8] = {{ 0.060376924683281, -0.00928875764374337, 0.00102893039315878, -9.55085533670919E-5, 7.58291285499256E-6, -5.01081057744427E-7, 2.4013441570345E-8, -5.1018669153887E-10 },
			      { 0.776823355931033, -0.119511285526388, 0.00974266885190267, -3.30447806384059E-4, -1.58051990661661E-5, 1.56022811158727E-6, 7.74602292865683E-8, -1.29646524960555E-8 },
			      { 6.66279971938553, -1.02504611065774, 0.0494758452357327, 2.44217481700129E-4, -7.32728109752881E-6, -2.507344770642E-6, -3.02786290067014E-7, -9.28536484109606E-9 }};
	constant__
	TYPE__ F_CP[10] = { 0.199999999999986, -0.142857142854412, 0.0555555554649585, -0.015151513983854, 0.00320512054753924, -5.55526624875562E-4, 8.16324851790106E-5, -1.03463270693454E-5, 1.09552870123182E-6, -7.6091148609885E-8 };
	// BOOST_GPU_ENABLED
	// static void evaluate(TYPE__ X, TYPE__ *R, TYPE__ *f) {
	//     evaluate1(X, R, R_CP);
	//     evaluate1(X, &f[1], F_CP);
	// }
	template<class A, typename T>
	BOOST_GPU_ENABLED
	static void evaluate(TYPE__ X, A &R, TYPE__ &f1,
			     const T &thread) {
	    polynomial::evaluate(R, X, R_CP, thread);
	    if (master(thread)) evaluate1(X, f1, F_CP);
	}
    }

    namespace x1to3 {
	constant__
	TYPE__  R_CP[3][10] = {{ 0.0452578254679079, -0.0061121934982509, 6.03194524398109E-4, -5.13430986707889E-5, 3.86118485517386E-6, -2.60122498274734E-7, 1.56592951656828E-8, -6.55098264095516E-10, 4.85300143926755E-12, 1.44687969563318E-12 },
			       { 0.573928229597613, -0.08487781203634, 0.00751549330892401, -3.89147469249594E-4, 9.923326947376E-7, 1.502366784525E-6, -6.745205954533E-8, -5.35281831445517E-9, 6.95964248788138E-10, 0.0 },
			       { 4.81234667357205, -0.824708946991557, 0.0504871005319119, 3.34761560498171E-5, -5.15964042227127E-5, -6.05865557561067E-6, -1.40971837780847E-7, 4.53631789436255E-8, 3.61058041895031E-9, -2.81496588401439E-10 }};
	constant__
	TYPE__ F_CP[12] = { 0.0529428148329709, -0.0323446977320647, 0.011384700111381, -0.00289955494844975, 5.83701488646511E-4, -9.74574633246452E-5, 1.39195169625425E-5, -1.73806555021045E-6, 1.92804632038796E-7, -1.92514145088973E-8, 1.78157031325097E-9, -1.4804423107214E-10 };
	template<class A, typename T>
	BOOST_GPU_ENABLED
	static void evaluate(TYPE__ X, A &R, TYPE__ &f1,
			     const T &thread) {
	    TYPE__ y = X - 2.0;
	    polynomial::evaluate(R, y, R_CP, thread);
	    if (master(thread)) evaluate1(y, f1, F_CP);
	}
    }

    namespace x3to5 {
	constant__
	TYPE__ R_CP[3][10] = {{ 0.0350898470729044, -0.00421006346373634, 3.70326938096287E-4, -2.87048671521677E-5, 2.026002142457E-6, -1.229940017368E-7, 7.649155832025E-9, -4.66622033006074E-10, 1.44265709189601E-11, 0.0 },
			      { 0.431180523260239, -0.0593483267268959, 0.00530601428208358, -3.33739312603632E-4, 1.11788717230514E-5, 5.15021914287057E-7, -7.95045680685193E-8, 2.15971131403034E-9, 1.97549041402552E-10, -2.65526039155651E-11 },
			      { 3.36412312243724, -0.624498381002855, 0.0489665802767005, -6.14120969315617E-4, -1.05296443527943E-4, -3.05512456576552E-6, 6.40574545989551E-7, 4.42413039572867E-8, -4.1642322978228E-9, -3.92833750584041E-10 }};
	constant__
	TYPE__ F_CP[12] = { 0.0175257821619922, -0.00866415899015349, 0.00264584212770942, -6.10614987696677E-4, 1.14504948737066E-4, -1.81353179793672E-5, 2.48749160874431E-6, -3.00873821471489E-7, 3.25336816562485E-8, -3.18111322308846E-9, 2.89147476459092E-10, -2.36788772599074E-11 };
	template<class A, typename T>
	BOOST_GPU_ENABLED
	static void evaluate(TYPE__ X, A &R, TYPE__ &f1, const T & thread) {
	    TYPE__ y = X - 4.0;
	    polynomial::evaluate(R, y, R_CP, thread);
	    if (master(thread)) evaluate1(y, f1, F_CP);
	}
    }

    namespace x5to10 {
	constant__
	TYPE__ R_CP[3][14] = {{ 0.0239112249488821, -0.00239864911618015, 1.77930381549953E-4, -1.08656008854077E-5, 7.34609900170759E-7, -4.65449552856856E-8, 1.39032379769474E-9, -4.23879635610964E-11, 1.31541892704E-11, -6.264613873998E-13, -6.736701449826E-14, 7.11884203790984E-16, 5.74429401360115E-16, 0.0 },
			      { 0.275969447451882, -0.0323845035189063, 0.00269443611274173, -1.68211091755327E-4, 9.68100917793911E-6, -4.40252529648056E-7, -1.06656985608482E-8, 2.73681574882729E-9, -2.492175211635E-11, -5.293620408757E-12, -8.595618132088E-13, 6.99375313934242E-15, 1.1346409620912E-14, 0.0 },
			      { 1.7358283175543, -0.32233505127086, 0.0351246831672571, -0.00183917335649633, -3.14232666170796E-5, 8.97224398620038E-6, 8.769581622041E-9, -6.205591993923E-8, 1.663165279876E-9, 3.917984522103E-10, -2.309293727603E-11, -1.985141104444E-12, 1.84955640200794E-13, 6.66339416996191E-15 }};
	constant__
	TYPE__ W0_CR[7] = { -3.1501078774085E-6, -0.49984072848436, 0.24645596956002, -0.32883030418398, 0.53689283271887, -0.69955602298985, 0.46897511375022 };
	template<class A, typename T>
	BOOST_GPU_ENABLED
	static void evaluate(TYPE__ X, A &R, A &W,
			     const TYPE__ e, const T & thread) {
	    TYPE__ y = X - 7.5;
	    polynomial::evaluate(R, y, R_CP, thread);
	    // }
	    if (master(thread)) W[0] = e*polynomial::inverse(W0_CR, X);
	}
    }

    namespace x10to15 {
	constant__
	TYPE__ R_CP[3][14] = {{ 0.0153435577063174, -0.00118587782909876, 8.04427643593792E-5, -4.14222202791434E-6, 1.29454702851936E-7, -3.55578027040563E-9, 7.11372337079797E-10, -8.67286219861085E-11, 1.882093066702E-12, 5.379885121517E-13, -4.084026087887E-14, -2.77189767070441E-15, 4.4213300128309E-16, 0.0 },
			      { 0.165077877454402, -0.0145361636398178, 0.00111888323089714, -6.63678914608681E-5, 2.51020389884249E-6, -5.59678576054633E-8, 7.08558457182751E-9, -1.13413908163831E-9, 4.798806828724E-11, 6.642452485783E-12, -8.579165965128E-13, -1.08257654410279E-14, 6.85146742119357E-15, 0.0 },
			      { 0.777203937334739, -0.10144965289945, 0.0113901881430697, -0.00101308723606946, 6.27436306915967E-5, -1.57938204115055E-6, -1.567725007761E-7, 2.064664199164E-8, -7.718080513708E-10, -5.625235879301E-11, 8.654129268056E-12, -3.157134329361E-13, -2.73458804864628E-14, 3.20622388697743E-15 }};
	constant__
	TYPE__ W0_CR[4] = { -2.1916512131607E-5, -0.49893752514047, 0.22991849164985, -0.18784686463512 };
	template<class A, typename T>
	BOOST_GPU_ENABLED
	static void evaluate(TYPE__ X, A &R, A &W,
			     const TYPE__ e, const T &thread) {
	    TYPE__ y = X - 12.5;
	    polynomial::evaluate(R, y, R_CP, thread);
	    if (master(thread)) W[0] = e*polynomial::inverse(W0_CR, X);
	}
    }

    namespace x15to20 {
	constant__
	TYPE__ R_CP[3][7] = {{ -2079.70687843258, 243.517435690398, -17.3209218219175, 0.781425144913975, -0.0234112415981143, 3.57901398988359E-4, -2.43270989903742E-6 },
			     { 33520.2872835409, -2366.59637247087, 107.037141010778, -3.0933761873188, 0.0349187925428138, -2.62627010965435E-4, 0.0 },
			     { -6881.45821789955, 404.996712650414, -18.4338896480695, -0.783503697918455, -0.0287029400759565, 9.31856404738601E-5, 0.0 }};
	constant__
	TYPE__ R_CR[3][3] = {{ 9824.41363463929, -19761.1541576986, 0.0 },
			     { -291532.335433779, 1411295.05262758, -2916691.1368102 },
			     { 51149.8390849158, -189829.509315154, 0.0 }};
	constant__
	TYPE__ W0_CR[3] = { -6.0156581186481E-5, -0.4969524146449, 0.1962326414943 };
	template<class A, typename T>
	BOOST_GPU_ENABLED
	static void evaluate(TYPE__ X, A &R, A &W,
			     const TYPE__ e, const T &thread) {
	    polynomial::evaluate(R, X, R_CP, thread);
	    polynomial::inverse0(R, R_CR, X, thread);
	    scale_add(e, R, X, R3, thread);
	    if (master(thread)) W[0] = e*polynomial::inverse(W0_CR, X);
	}
    }

    namespace x20to33 {
	constant__
	TYPE__ R_CP[3][5] = {{ 164.931462413877, -18.8336409225481, 1.31099142238996, -0.0500929599665316, -4.97561537069643E-4 },
			     { 1522.31757709236, -165.426392885291, 11.3691058739678, -0.517373211334924, -0.00448218898474906 },
			     { 2698.31813951849, -357.615122086961, 17.3639054044562, -1.77293428863008, -0.0138368602394293 }};
	constant__
	TYPE__ R_CR[3][1] = {{ -660.344754467191 },
			     { -6309.09125686731 },
			     { -14573.4701095912 }};
	constant__
	TYPE__ W0_CR[3] = { -6.0156581186481E-5, -0.4969524146449, 0.1962326414943 };
	template<class A, typename T>
	BOOST_GPU_ENABLED
	static void evaluate(TYPE__ X, A &R, A &W,
			     const TYPE__ e, const T &thread) {
	    polynomial::evaluate(R, X, R_CP, thread);
	    polynomial::inverse0(R, R_CR, X, thread);
	    scale_add(e, R, X, R3, thread);
	    if (master(thread)) W[0] = e*polynomial::inverse(W0_CR, X);
	}
    }

    namespace x33to47 {
	constant__
	TYPE__  R_CP[3][3] = {{ -3994.33696473658, 321.318352526305, -7.39058467995275 },
			      { -38686.2867311321, 3135.69966333873, -73.8726243906513 },
			      { -128094.577915394, 10441.2168692352, -263.750565461336 }};
	constant__
	TYPE__ W_CP[3][4] = {{ 0.0, 0.0, 0.0, 0.0 },
			     { 38079.4303087338, -2919.80647450269, 61.5072615497811, 0.0 },
			     { -1677.87926005344, 192.977367967984, -8.30661900042651, 0.152258947224714 }};
	template< class A, typename T>
	BOOST_GPU_ENABLED
	static void evaluate(TYPE__ X, A &R, A &W,
			     const TYPE__ e, const TYPE__ w0,
			     const T &thread) {
	    polynomial::evaluate2(X, R, W, R_CP, W_CP, thread);
	    scale_add(e, R, W, X, R3, w0, W3, thread);
	    if (master(thread)) W[0] -= (W[1] + W[2]);
	}
    }

#define F(i) W[(i)+1]

    BOOST_GPU_ENABLED
    static void weights3(TYPE__ (&R)[3], TYPE__ (&W)[3], const serial_tag&) {
	TYPE__ a[2];

	a[1] = F(1) - R[0]*F(0);
	a[0] = F(0) - R[0]*W[0];

	W[1] = (R[2]*a[0] - a[1])/((R[2] - R[1])*(R[1] - R[0]));
	W[2] = (a[1] - R[1]*a[0])/((R[2] - R[1])*(R[2] - R[0]));

	W[0] -= (W[1] + W[2]);
    }

    template<class A, typename T>
    BOOST_GPU_ENABLED
    static void weights3(A &R, A &W, const T &thread) {
	if (index(thread) >= 2) return;
	TYPE__ a;

	a = F(thread) - R[0]*((thread == 0) ? W[0] : F(0));
	F(thread) = a;
	a = (R[!(thread%2) + 1]*F(0) - F(1));
	a *= (1 - 2*TYPE__(thread));

	W[thread+1] = a/((R[2] - R[1])*(R[thread+1] - R[0]));

	if (master(thread)) W[0] -= (W[1] + W[2]);
    }

    template<class A, typename T>
    BOOST_GPU_ENABLED
    static void evaluate_(TYPE__ X, A &R, A &W, const T &thread) {

	if (index(thread) >= N) return;

	const TYPE__ e = exp(-X);
	const TYPE__ W0 = sqrt(PIE4/X); 

	if (X > 47.0) {
	    divide(R, X, R3, thread);
	    scale(W, W0, W3, thread);
	    if (master(thread)) W[0] -= (W[1] + W[2]);
	    return;
	}
	else if (X > 33.0) {
	    x33to47::evaluate(X, R, W, e, W0, thread);
	}
	else if (X > 5.0) {
	    if (X > 20.0) x20to33::evaluate(X, R, W, e, thread);
	    else if (X > 15.0) x15to20::evaluate(X, R, W, e, thread);
	    else if (X > 10.0) x10to15::evaluate(X, R, W, e, thread);
	    else x5to10::evaluate(X, R, W, e, thread);
	    if (master(thread)) {
		W[0] += W0;
		F(0) = (W[0] - e)/(X + X);
		F(1) = (3*F(0) - e)/(X + X);
	    }
	}
	else if (X > 3e-7) {
	    if (X > 3.0) x3to5::evaluate(X, R, F(1), thread);
	    else if (X > 1.0) x1to3::evaluate(X, R, F(1), thread);
	    else x0to1::evaluate(X, R, F(1), thread);
	    if (master(thread)) {
		F(0) = ((X + X)*F(1) + e)/3.0;
		W[0] = (X + X)*F(0) + e;
	    }
	}
	else {
	    x0::evaluate(X, R, W, thread);
	}

	change_variable<N>(R, thread);
	if (X > 33.0) return;
	if (X > 3e-7) weights3(R, W, thread);

#undef F
	
    }


    BOOST_GPU_ENABLED
    inline void evaluate(TYPE__ X, TYPE__ (&R)[3], TYPE__ (&W)[3]) {
	evaluate_(X, R, W, serial_tag());
    }

    template<typename T>
    BOOST_GPU_ENABLED
    inline void evaluate(TYPE__ X, TYPE__ *R, TYPE__ *W, const T &thread) {
	evaluate_(X, R, W, thread);
    }

} // namespace roots3
