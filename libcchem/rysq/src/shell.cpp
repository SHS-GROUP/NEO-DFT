#include <algorithm>

#include "rysq-core.hpp"

using namespace rysq;

// Shell::Shell(rysq::type type, int K, const double *a, int nc, const double* C[])
//     :    shell(type, K)
// }

// Shell::Shell(int L,  int K, const double *a, const double *C) {
//     const double *C_[] = { C };
//     initialize(rysq::type(L), K, a, 1, C_);
//     //impl_ = new Shell::Impl(rysq::type(L), K, a, 1, C_);
// }

// Shell::Shell(int K, const double *a, const double *Cs, const double *Cp) {
//     const double *C_[] = { Cs, Cp };
//     //impl_ = new Shell::Impl(rysq::SP, K, a, 2, C_);
//     initialize(rysq::SP, K, a, 2, C_);
// }

static char get_symbol(int L) {
    static const char symbols[] = "spd";
    return (L < 3) ? symbols[L] : char('f' + L - 3);
}

static void set_symbol(rysq::type type, char *symbol) {
    if (type < 0) {
	for (int i = 0; i <= std::abs(type); ++i) {
	    symbol[i] = get_symbol(i);
	}
	symbol[std::abs(type)+1] = '\0';
    }
    else {
	symbol[0] = get_symbol(type);
	symbol[1] = '\0';
    }
}


void Shell::initialize(const double *a, const double* const C[]) {

    this->data_ = new double[K*(1 + nc)];
    this->a_ = this->data_;
    std::copy(a, a + K, this->a_);

    this->C_[0] = this->C_[1] = NULL; 
    for (int i = 0; i < nc; ++i) {
	this->C_[i] = this->data_ + (i+1)*K;
	std::copy(C[i], C[i] + K, this->C_[i]);
    }
    
    set_symbol(this->type, this->symbol_);

}

