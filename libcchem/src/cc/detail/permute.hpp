#ifndef CC_DETAIL_PERMUTE_HPP
#define CC_DETAIL_PERMUTE_HPP

#include "cc/cc.hpp"

#include "array/permute.hpp"
#include <boost/multi_array.hpp>

namespace cchem {
namespace cc {
namespace detail {

    template<int O>
    void permute(const Array &A, Array &T);

    template<>
    void permute<0321>(const Array &A, Array &T) {
	assert(A.N == 4);
	assert(T.N == 4);

	int ni = A.shape()[0];
	int nj = A.shape()[1];
	int nk = A.shape()[2];
	int nl = A.shape()[3];

	boost::multi_array<double,3> a(boost::extents[nl][nj][ni]);
	boost::multi_array<double,3> t(boost::extents[nj][nl][ni]);

	for (int k = 0; k < nk; ++k) {
	    size_t r1[] = { 0, 0, k, 0 };
	    size_t r2[] = { ni, nj, k+1, nl };
	    A.get(a.data(), r1, r2);
	    for (int l = 0; l < nl; ++l) {
		for (int j = 0; j < nj; ++j) {
		    for (int i = 0; i < ni; ++i) {
			t[j][l][i] = a[l][j][i];
		    }
		}
	    }
	    std::swap(r1[1], r1[3]);
	    std::swap(r2[1], r2[3]);
	    // std::cout << r2[0] << " " << r2[1] << " "
	    // 	      << r2[2] << " " << r2[3] << std::endl;
	    T.put(t.data(), r1, r2);
	}
    }

    template<>
    void permute<2103>(const Array &A, Array &T) {
	assert(A.N == 4);
	assert(T.N == 4);

	int ni = A.shape()[0];
	int nj = A.shape()[1];
	int nk = A.shape()[2];
	int nl = A.shape()[3];

	boost::multi_array<double,3> a(boost::extents[nk][nj][ni]);
	boost::multi_array<double,3> t(boost::extents[ni][nj][nk]);

	for (int l = 0; l < nl; ++l) {
	    size_t r1[] = { 0, 0, 0, l };
	    size_t r2[] = { ni, nj, nk, l+1 };
	    A.get(a.data(), r1, r2);
	    for (int k = 0; k < nk; ++k) {
		for (int j = 0; j < nj; ++j) {
		    for (int i = 0; i < ni; ++i) {
			t[i][j][k] = a[k][j][i];
		    }
		}
	    }
	    std::swap(r1[2], r1[0]);
	    std::swap(r2[2], r2[0]);
	    T.put(t.data(), r1, r2);
	}
    }

    template<>
    void permute<0213>(const Array &A, Array &T) {
	assert(A.N == 4);
	assert(T.N == 4);

	int ni = A.shape()[0];
	int nj = A.shape()[1];
	int nk = A.shape()[2];
	int nl = A.shape()[3];

	boost::multi_array<double,3> a(boost::extents[nk][nj][ni]);
	boost::multi_array<double,3> t(boost::extents[nj][nk][ni]);

	for (int l = 0; l < nl; ++l) {
	    size_t r1[] = { 0, 0, 0, l };
	    size_t r2[] = { ni, nj, nk, l+1 };
	    A.get(a.data(), r1, r2);
	    for (int k = 0; k < nk; ++k) {
		for (int j = 0; j < nj; ++j) {
		    for (int i = 0; i < ni; ++i) {
			t[j][k][i] = a[k][j][i];
		    }
		}
	    }
	    std::swap(r1[2], r1[1]);
	    std::swap(r2[2], r2[1]);
	    T.put(t.data(), r1, r2);
	}
    }

} // namespace detail
} // namespace cc
} // namespace cchem

#endif /* CC_DETAIL_PERMUTE_HPP */
