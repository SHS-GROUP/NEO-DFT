#ifndef INTEGRAL_SCREENING_HPP
#define INTEGRAL_SCREENING_HPP

#include "core/integral.hpp"

namespace integral {

struct screening {
    screening() : cutoff_(0) {}
    screening(const Matrix &K, double cutoff = 1e-10)
	: K_(K), cutoff_(cutoff) {}
    bool test(double value) const {return (std::abs(value) > cutoff_); }
    bool test(int i, int j, int k, int l) const {
	if (K_.size1() <= 0) return true;
	return (K_(i,j)*K_(k,l) > cutoff_);
    }
private:
    const Matrix K_;
    double cutoff_;
};

}

#endif // INTEGRAL_SCREENING_HPP
