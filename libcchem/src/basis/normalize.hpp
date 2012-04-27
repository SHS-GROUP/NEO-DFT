#ifndef BASIS_NORMALIZE_HPP
#define BASIS_NORMALIZE_HPP

#include "basis/basis.hpp"
#include "foreach.hpp"

namespace basis {

    template<class V>
    void normalize_vector(const Basis &basis, V &v) {
	typename V::iterator it = v.begin();
	foreach (const Shell &shell, basis) {
	    foreach (const Shell::Function &f, shell.functions()) {
		*(it++) *= f.C;
	    }
	}
    }

    template<class A>
    void normalize_matrix(const Basis &basis, A &a) {
	std::vector<double> c(basis.size(), 1);
	normalize_vector(basis, c);
	for (size_t j = 0; j < c.size(); ++j) {
	    for (size_t i = 0; i < c.size(); ++i) {
		a(i,j) *= c[i]*c[j];
	    }
	}
    }

    template<class A>
    void normalize_matrix1(const Basis &basis, A &m) {
	std::vector<double> c(basis.size(), 1);
	normalize_vector(basis, c);
	for (size_t j = 0; j < size_t(m.size2()); ++j) {
	    for (size_t i = 0; i < c.size(); ++i) {
		// std::cout << c[i] << std::endl;
		m(i,j) *= c[i];
	    }
	}
    }

}

#endif // BASIS_NORMALIZE_HPP
