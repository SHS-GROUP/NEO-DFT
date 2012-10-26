#include "basis/shell.hpp"

#include "utility/numeric.hpp"
#include <math.h>
#include <vector>
#include <algorithm>
#include <iostream>

double basis::shell_data::Function::coefficient(int i, int j, int k) {
    using utility::factorial2;
    int n = std::max(i + j + k, 1);
    int nx = std::max(i, 1);
    int ny = std::max(j, 1);
    int nz = std::max(k, 1);
    double denom = factorial2(2*nx -1 )*factorial2(2*ny -1 )*factorial2(2*nz -1 );
    return sqrt(factorial2(2*n - 1)/denom);
}

std::vector<basis::shell_data::Function> generate_functions(int L) {
    std::vector<basis::shell_data::Function> result;
    for (int i = L; i >= 0; --i) {
	for (int k = 0; k <= (L-i)/2; ++k) {
	    int j = L - i - k;
	    if (i < j || j < k) continue;
	    utility::permute(i, j, k, result);
	}
    }
    return result;
}

/** @see basis.h */
basis::shell_data::shell_data(int K, double *exps, double *C, int L)
    : type(shell::type(L))
{
    K_ = K;
    a_.assign(exps, exps + K);

    nc_ = 1;
    C_.assign(1, std::vector<double>(K));
    C_[0].assign(C, C + K);

    L_ = L;
    Lmin_ = L;

    size_ = (((L+1)*(L+1) + (L+1))/2);

    functions_ = generate_functions(L);
}

/** @see basis.h */
basis::shell_data::shell_data(int K, double *exps, double *Cs, double *Cp)
    : type(shell::type(-1))
{
    K_ = K;

    a_.assign(exps, exps+ K);

    nc_ = 2;
    C_.assign(2, std::vector<double>(K));
    C_[0].assign(Cs, Cs + K);
    C_[1].assign(Cp, Cp + K);

    L_ = 1;
    Lmin_ = 0;
    size_ = 4;

    for (int l = Lmin_; l <= L_; ++l) {
	std::vector<basis::shell_data::Function> list = generate_functions(l);
	functions_.insert(functions_.end(), list.begin(), list.end());
    }

}

bool basis::shell_data::operator==(const basis::shell_data &shell) const {
    if (K_ != shell.K_) return false;
    else if (type != shell.type) return false;
    else if (nc_ != shell.nc_) return false;
    else if (a_ != shell.a_) return false;
    else if (C_ != shell.C_) return false;

    return true;
}

std::ostream& basis::operator<<(std::ostream& os, const Shell &shell) {
    int type = shell.data()->type;
    static const char symbol[] = "spd";
    if (type < 0)
	return (os << std::string(symbol, 1 + -type));
    if (type < 3)
	return (os << symbol[type]);
    return (os << char('f'+type-3));
}

std::ostream& basis::operator<<(std::ostream &output, const basis::shell_data &s) {
    output << "ptr = " << &s << "\tK = " << s.K() << "\t"
	   << "L = " << s.type << std::endl;

//     for (unsigned int i = 0; i < s.functions.size(); ++i) {
// 	output << "("
// 	       << (s.functions.at(i))(0) << ", "
// 	       << (s.functions.at(i))(1) << ", "
// 	       << (s.functions.at(i))(2) << ")"
// 	       << "\t" << s.functions.at(i).C
// 	       << std::endl;
//     }

    for (int i = 0; i < s.K(); ++i) {
	output.width(13);
        output.precision(8);
	output << s(i);

	for (int j = 0; j < s.nc(); ++j) {
	    output.width(13);
	    output.precision(8);
	    output << '\t' << s(i, j);
	}

	output << std::endl;
    }

    return output;
}

