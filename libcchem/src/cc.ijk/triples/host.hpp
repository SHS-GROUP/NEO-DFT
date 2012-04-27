#ifndef CC_TRIPLES_HOST_HPP
#define CC_TRIPLES_HOST_HPP

#include "cc/triples/triples.hpp"
#include "cc/triples/thread.hpp"

#include <cuda.h>
#include "blas.hpp"

namespace cc {
namespace triples {
namespace detail {


    struct Host
	: detail::Thread< array::adapter<double,3> >
    {
	typedef array::adapter<double,3> T3;
	typedef detail::Thread< array::adapter<double,3> > Thread;
    private:
	size_t num_threads_;
	Matrix vv_[2];

    public:
	template<size_t N, typename U>
	static array::adapter<double,N>
	make_array(const U &shape, double *data = 0) {
	    size_t size = 1;
	    for (size_t i = 0; i < N; ++i)
		size *= shape[i];
	    if (!data)
		data = allocate<double>(size);
	    return array::adapter<double,N>(data, shape);
	}
	template<typename T>
	static T* allocate(size_t size) {
	    return new T[size];
	}
	template<typename T>
	static void free(T *ptr) {
	    delete ptr;
	}

    public:
	Host(const Thread &thread)
	    : Thread(thread)
	{
	    //size_t no = Thread::no_;
	    size_t nv = Thread::nv_;
	    vv_[0].resize(nv, nv, false);
	    vv_[1].resize(nv, nv, false);
	}

	void ijk1_impl(const Data &data,
		       const Matrix &tab, const Matrix &taib,
		       const Matrix &Vjk, const Matrix &Vkj,
		       T3 &t3,
		       int i, int j, int k,
		       std::pair<bool,int> terms,
		       Queue &queue) {

	    namespace ublas = boost::numeric::ublas;
	    Layout layout;

	    int no = data.no;
	    int nv = data.nv;

	    while (true) {
		BOOST_AUTO(task, queue.advance());
		if (task.size < 1) break;
		int c = task.begin;
		int n = task.size;

		double b = 0;
		Matrix &t3_c = vv_[0];
		Matrix &t3_b = vv_[1];

		if (terms.first) {
		    Matrix &Vab = vv_[1];
		    data["V(abci)"]
		    	.get(Vab.data().begin(),
			     index<4>(0,0,c,k), index<4>(nv,nv,c+n,k+1));
		    blas::gemm(1, tab, Vab, b, t3_c);
		    blas::gemm(1, Vab, tab, 1, t3_c);
		    b = 1;
		}

		BOOST_AUTO(tai, ublas::project(taib,
					       ublas::range(0,nv),
					       ublas::range(c*no, (c+1)*no)));
		if (terms.second == 2) {
		    BOOST_AUTO(tai, ublas::project(taib,
		    				   ublas::range(0,nv),
		    				   ublas::range(c*no, (c+1)*no)));
		    blas::gemm(-1, tai, Vjk, b, t3_c);
		}

		if (terms.second != 1) {
		    BOOST_AUTO(t3_,
			       ublas::make_matrix(nv, nv, t3[c].origin(), layout));
		    ublas::noalias(t3_) += t3_c;
		}
		// t3(*,c,*) 
		else {
		    blas::gemm(-1, tai, Vkj, 0, t3_b);
		    for (int b = 0; b < nv; ++b) {
			BOOST_AUTO(t3_bc,
				   ublas::make_vector(nv, t3[b][c].origin()));
			ublas::noalias(t3_bc) += ublas::column(t3_b,b);
		    }
		}

	    } // while

	}

    };

}
}
}

#endif // CC_TRIPLES_HOST_HPP

