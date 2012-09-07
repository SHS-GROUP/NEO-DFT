#ifndef ARRAY_PERMUTE_HPP 
#define ARRAY_PERMUTE_HPP

#include <boost/mpl/vector_c.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/utility/swap.hpp>
#include <boost/typeof/typeof.hpp>

#include "array/adapter.hpp"
#include "array/detail/tile.hpp"
#include "array/detail/grid.hpp"

#ifdef __CUDACC__
#include "array/detail/device.hpp"
#endif // __CUDACC__

namespace array {
namespace detail {


    struct permute {

	template<size_t N, class A, class B, class G, class TB>
	BOOST_GPU_ENABLED
	void operator()(const A &a, B &b, boost::mpl::vector_c<int,1,0>,
			const G &g, const TB &tb, size_t ni, size_t nj) {
	    BOOST_STATIC_ASSERT(A::dimensionality == 2);
	    BOOST_STATIC_ASSERT(B::dimensionality == 2);
	    const int Ni = ni/N + (ni%N > 0);
	    const int Nj = nj/N + (nj%N > 0);

	    for (int j = g.y(); j < Nj; j += g.ny()) {
		range rj(j*N, min<int>((j+1)*N,nj));
		for (int i = g.x(); i < Ni; i += g.nx()) {
		    range ri(i*N, min<int>((i+1)*N,ni));
		    tile<typename B::element,N> bij;
		    bij.load(ri, rj, a, tb);
		    bij.transpose(rj, ri, b, tb);
		}
	    }
	}


	template<size_t B, class A, class G, class TB>
	BOOST_GPU_ENABLED
	void operator()(A &a, boost::mpl::vector_c<int,1,0>,
			const G &g, const TB &b, size_t n) {
	    BOOST_STATIC_ASSERT(A::dimensionality == 2);
	    typedef typename A::element T;
	    const int nB = n/B + (n%B > 0);

	    for (int j = g.y(); j < nB; j += g.ny()) {
		range rj(j*B, min<int>((j+1)*B,n));
		for (int i = g.x(); i <= j; i += g.nx()) {
		    range ri(i*B, min<int>((i+1)*B,n));
		    tile<T,B> aij, aji;
		    aij.load(ri, rj, a, b);
		    aji.load(rj, ri, a, b);
		    aij.transpose(rj, ri, a, b);
		    if (i == j) continue;
		    aji.transpose(ri, rj, a, b);
		}
	    }
	}

	template<class A, class B, class G, class T>
	BOOST_GPU_ENABLED
	void operator()(const A &a, B &b, boost::mpl::vector_c<int,1,0>,
			const G &g, const T &tb) {
	    BOOST_STATIC_ASSERT(A::dimensionality == 2);
	    BOOST_STATIC_ASSERT(B::dimensionality == 2);
	    BOOST_ASSERT(a.shape()[0] == b.shape()[1]);
	    BOOST_ASSERT(b.shape()[0] == a.shape()[1]);
	    operator()<16>(a, b, boost::mpl::vector_c<int,1,0>(), g, tb,
			   a.shape()[0], a.shape()[1]);
	    // for (size_t j = 0; j < a.shape()[1]; ++j) {
	    //     for (size_t i = 0; i < j; ++i) {
	    // 	std::swap(a[i][j], a[j][i]);
	    //     }
	    // }
	}

	template<class A, class G, class T>
	BOOST_GPU_ENABLED
	void operator()(A &a, boost::mpl::vector_c<int,1,0>,
			const G &g, const T &b) {
	    BOOST_STATIC_ASSERT(A::dimensionality == 2);
	    BOOST_ASSERT(a.shape()[0] == a.shape()[1]);
	    operator()<16>(a, boost::mpl::vector_c<int,1,0>(), g, b, a.shape()[0]);
	    // for (size_t j = 0; j < a.shape()[1]; ++j) {
	    //     for (size_t i = 0; i < j; ++i) {
	    // 	std::swap(a[i][j], a[j][i]);
	    //     }
	    // }
	}

	template<class A, class G, class B>
	BOOST_GPU_ENABLED
	void operator()(A &a, boost::mpl::vector_c<int,1,0,2>,
			const G &g, const B &b) {
	    BOOST_STATIC_ASSERT((A::dimensionality == 3));
	    BOOST_ASSERT(a.shape()[1] == a.shape()[2]);
	    for (size_t k = g.z(); k < a.shape()[0]; k += g.nz()) {
		BOOST_AUTO(ak, a[k]);
		operator()(ak, boost::mpl::vector_c<int,1,0>(), g, b);
	    }
	}

	template<class A, class G, class B>
	BOOST_GPU_ENABLED
	void operator()(A &a, boost::mpl::vector_c<int,0,2,1>,
			const G &g, const B &b) {
	    BOOST_STATIC_ASSERT((A::dimensionality == 3));
	    // std::cout << a.shape()[1] << " " << a.shape()[2]
	    // 	      << " " << a.shape()[0] << std::endl;
	    BOOST_ASSERT(a.shape()[0] == a.shape()[1]);
	    for (size_t j = g.y(); j < a.shape()[0]; j += g.ny()) {
		for (size_t i = g.x(); i < j; i += g.nx()) {
		    detail::swap(a[i][j].begin(), a[i][j].end(), a[j][i].begin(), b);
		}
	    }
	}

	// template<class A, class G, class B>
	// void operator()(A &a, boost::mpl::vector_c<int,2,1,0>,
	// 		const G &g, const B &b) {
	//     BOOST_ASSERT(a.shape()[0] == a.shape()[2]);
	//     size_t K = a.shape()[1];
	//     for (size_t j = g.y(); j < a.shape()[2]; j += g.ny()) {
	// 	for (size_t k = 0; k < K; ++k) {
	// 	    for (size_t i = g.x(); i < j; i += g.nx()) {
	// 		detail::swap(a[j][k][i], a[i][k][j]);
	// 	    }
	// 	}
	//     }
	// }

    };

}
}

namespace array {

    template<int I, int J, class A, class B>
    void permute(const A &a, B &b) {
	BOOST_STATIC_ASSERT((detail::rank<A>::value == 2 &&
			     detail::rank<B>::value == 2));
	using detail::grid;
	adapter<const typename A::element,2> a_(a.origin(), a.shape());
	adapter<typename B::element,2> b_(b.origin(), b.shape());
	detail::permute()(a_, b_, boost::mpl::vector_c<int,I,J>(),
			 grid::serial(), grid::serial());
    }

    template<int I, int J, class A>
    void permute(A &a) {
	BOOST_STATIC_ASSERT((detail::rank<A>::value == 2));
	using detail::grid;
	adapter<typename A::element,2> a_(a.origin(), a.shape());
	detail::permute()(a_, boost::mpl::vector_c<int,I,J>(),
			  grid::serial(), grid::serial());
    }

    template<int I, int J, int K, class A>
    void permute(A &a, int id, size_t count) {
	BOOST_STATIC_ASSERT((detail::rank<A>::value == 3));
	using detail::grid;
	adapter<typename A::element,3> a_(a.origin(), a.shape());
	typename boost::mpl::if_c<(K == 2),
	    grid::Z, grid::Y>::type g(id, count);
	detail::permute()(a_, boost::mpl::vector_c<int,I,J,K>(),
			  g, grid::serial());
    }

    template<int I, int J, int K, class A>
    void permute(A &a) {
	permute<I,J,K>(a, 0, 1);
    }

#ifdef __CUDACC__
    template<int I, int J, int K, typename T, size_t N>
    void permute(adapter<T, N, device_tag> &a, int id, size_t count) {
	if (id != 0) return;
	//std::cout << "permute" << std::endl;
	//BOOST_STATIC_ASSERT((0));

	dim3 grid((a.shape()[1]+15)/16, (a.shape()[2]+15)/16);
	dim3 block(16,16);

	if (I == 0 && J == 2 && K == 1) {
	    grid = dim3(a.shape()[1], a.shape()[2]);
	    block = dim3(64);
	}

	(detail::device::kernel(grid, block)
	 (detail::permute(), a, boost::mpl::vector_c<int,I,J,K>()));
    }
#endif // __CUDACC__


}

#endif // ARRAY_PERMUTE_HPP
