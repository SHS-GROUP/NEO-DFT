#ifndef BASIS_NORMALIZE_HPP
#define BASIS_NORMALIZE_HPP

#include "basis/basis.hpp"
#include "foreach.hpp"

#include <assert.h>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>

namespace basis {

    // template<class V>
    // void normalize_vector(const Basis &basis, V &v) {
    // 	typename V::iterator it = v.begin();
    // 	foreach (const Shell &shell, basis) {
    // 	    foreach (const Shell::Function &f, shell.functions()) {
    // 		*(it++) *= f.C;
    // 	    }
    // 	}
    // }

    template<int K, typename T, class A>
    static void normalize(const std::vector<T> &N, A &a) {
	BOOST_STATIC_ASSERT((0 <= K && K <= 2));
	if (K == 0 || K == 2) assert(a.size1() == N.size());
	if (K == 1 || K == 2) assert(a.size2() == N.size());
	for (size_t j = 0; j < a.size2(); ++j) {
	    double Nj = ((K == 1 || K == 2) ? N[j] : 1);
	    for (size_t i = 0; i < a.size1(); ++i) {
		double Ni = ((K == 0 || K == 2) ? N[i] : 1);
		a(i,j) *= (Nj*Ni);
	    }
	}
    }

    // template<typename T, class A>
    // static void normalize(const std::vector<T> &N,
    // 			  boost::numeric::ublas::matrix_expression<A> &a) {
    // 	assert(a().size1() != N.size());
    // 	for (size_t j = 0; j < a().size2(); ++j) {
    // 	    BOOST_AUTO(aj, column(a(),j));
    // 	    for (size_t i = 0; i < a().size1(); ++i) {
    // 		aj[i] *= N[i];
    // 	    }
    // 	}
    // }

    // template<class A>
    // void normalize_matrix(const Basis &basis, A &a) {
    // 	std::vector<double> c(basis.size(), 1);
    // 	normalize_vector(basis, c);
    // 	for (size_t j = 0; j < c.size(); ++j) {
    // 	    for (size_t i = 0; i < c.size(); ++i) {
    // 		a(i,j) *= c[i]*c[j];
    // 	    }
    // 	}
    // }

    // template<class A>
    // void normalize_matrix1(const Basis &basis, A &m) {
    // 	std::vector<double> c(basis.size(), 1);
    // 	normalize_vector(basis, c);
    // 	for (size_t j = 0; j < size_t(m.size2()); ++j) {
    // 	    for (size_t i = 0; i < c.size(); ++i) {
    // 		// std::cout << c[i] << std::endl;
    // 		m(i,j) *= c[i];
    // 	    }
    // 	}
    // }

}

#endif // BASIS_NORMALIZE_HPP
