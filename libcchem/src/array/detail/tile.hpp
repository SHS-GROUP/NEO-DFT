#ifndef ARRAY_DETAIL_TILE_HPP 
#define ARRAY_DETAIL_TILE_HPP

#include "array/detail/utility.hpp"
#include <boost/config.hpp>
#include <boost/typeof/typeof.hpp>

namespace array {
namespace detail {

    /**
       Array tile
       @tparam T value type
       @tparam N tile size
    */
    template<typename T, size_t N>
    struct tile {

	/**
	   Load tile from a[rj][ri].
	   Thread block argument allows to map I/O to 2D block.
	   @param ri range of i
	   @param rj range of j
	   @param a source array, a[rj][ri]
	   @param b thread block
	*/
	template<class R, class A, class B>
	BOOST_GPU_ENABLED
	void load(const R &ri, const R &rj, const A &a, const B &b) {
	    for (int j = b.y(); j < rj.size(); j += b.ny()) {
		typename A::const_reference aj =  a[j+rj.start()];
		for (int i = b.x(); i < ri.size(); i += b.nx()) {
		    data_[j][i] = aj[i+ri.start()];
		}
	    }
	}

	/**
	   transpose and store tile in a[rj][ri]
	*/
	template<class R, class A, class B>
	BOOST_GPU_ENABLED
	void transpose(const R &ri, const R &rj, A &a, const B &b) const {
	    for (int j = 0; j < rj.size(); ++j) {
		typename A::reference aj = a[j+rj.start()];
		for (int i = 0; i < ri.size(); ++i) {
		    aj[i+ri.start()] = data_[i][j];
		}
	    }
	}

	template<class R, class A>
	BOOST_GPU_ENABLED
	void load(const R &ri, const R &rj, const A &a) {
	    load(ri, rj, a, serial());
	}

	template<class R, class A>
	BOOST_GPU_ENABLED
	void transpose(const R &ri, const R &rj, A &a) const {
	    transpose(ri, rj, a, serial());
	}

    private:
	T data_[N][N];
    };

}
}

#endif // ARRAY_DETAIL_TILE_HPP
