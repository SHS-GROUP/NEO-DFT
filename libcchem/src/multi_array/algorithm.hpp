#ifndef MULTI_ARRAY_ALGORITHM_HPP
#define MULTI_ARRAY_ALGORITHM_HPP

#include <algorithm>
#include <boost/bind.hpp>
#include <boost/typeof/typeof.hpp>

namespace multi_array {
namespace algorithm {

namespace detail {

    template<class, class, class, class = void>
    struct functor;
	

    template<int Order>
    void permute(const size_t &i, size_t &j, size_t &k, size_t &l);

    template<>
    BOOST_GPU_ENABLED
    inline void permute<0231>(const size_t &i, size_t &j, size_t &k, size_t &l) {
	size_t t = j;
	j = k;
	k = l;
	l = t;
    }
	

    template<>
    BOOST_GPU_ENABLED
    inline void permute<0123>(const size_t &i, size_t &j, size_t &k, size_t &l) {
    }

    struct sequential {
	static const int x = 0, y = 0;
	static const int xsize = 1, ysize = 1;
    };

    template<int Order, class T0, class T1, class T>
    struct kernel2 {
	typedef functor<T, T0, T1> F;
	kernel2(T0 A, T1 B, F f = F())
	    : f(f), A(A), B(B) {}
	void operator()() {
	    operator()(sequential(), sequential());
	}
	template<class Grid, class Block>
	BOOST_GPU_ENABLED
	void operator()(const Grid &grid, const Block &block) {
	    for (size_t s = grid.x; s < size(0); s += grid.xsize) {
		for (size_t r = grid.y; r < size(1); r += grid.ysize) {
		    for (size_t q = 0; q < A.size[2]; ++q) {
			size_t b = q, c = r, d = s;
			permute<Order>(0,b,c,d);
			f(A[s][r][q].begin(), A[s][r][q].end(),
			  B[d][c][b].begin(), block);
		    }
		}
	    }
	}
	BOOST_GPU_ENABLED
	size_t size(int i) const { return A.size[i]; }
    private:
	F f; T0 A; T1 B;
    }; 

    template<class = void>
    struct copy;

    template<class T, class U>
    struct functor<copy<>, T, U> {
	template<class I, class O>
	void operator()(I begin, I end, O output, sequential) const {
	    std::copy(begin, end, output);
	}
    };

} // namespace detail
    
    template<int Order, typename T0, typename T1>
    detail::kernel2<Order, T0, T1, detail::copy<> >
    copy(T0 t0, T1 t1) {
	// typedef detail::copy<T0, T1> F;
	return detail::kernel2<Order, T0, T1, detail::copy<> >(t0, t1);
    }
	
} // namespace algorithm
} // namespace multi_array

#endif // MULTI_ARRAY_ALGORITHM_HPP
