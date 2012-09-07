#ifndef DFT_XC_KERNEL_HPP
#define DFT_XC_KERNEL_HPP

#include <math.h>
const double PI = 3.14159265358979323846;


struct Point {
    struct Density {
	double w, a, b, aa, bb, ab, dx, dy, dz;
    };
    struct XC2 {
	struct tuple { double V, dx, dy, dz, sigk, dk; };
	double Xa;
	tuple a, b;
    };
    struct XC {
	double Xa, Xg, Ec;
	double V, dx, dy, dz, sigk, dk;
    };
};


template<int N>
struct needs_gradient {
    struct gradient {
	static const int value = N;
    };
};

#endif // DFT_XC_KERNEL_HPP

