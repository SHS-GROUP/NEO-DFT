#ifndef EIGEN_SORT_HPP
#define EIGEN_SORT_HPP

#include <algorithm>
#include <boost/typeof/typeof.hpp>

namespace eigen {

    template<class V0, class A0>
    void sort(V0 &d, A0 &V) {
	int N = d.size();
	for (int i = 0; i < N-1; ++i) {
	    int k = i;
	    {
		BOOST_AUTO(p, d(i));
		for (int j = i+1; j < N; ++j) {
		    if (d(j) < p) { p = d(j); k = j; }
		}
	    }
	    if (k != i) {
		std::swap(d(k),d(i));
		for (int j = 0; j < N; ++j) std::swap(V(j,i), V(j,k));
	    }
	}
    }

}

#endif // EIGEN_SORT_HPP
