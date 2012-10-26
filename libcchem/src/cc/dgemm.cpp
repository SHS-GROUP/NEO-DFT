#include "cc/utility.hpp"
#include "cc/tensor.hpp"
#include "utility/timer.hpp"
#include "foreach.hpp"

#include "omp.hpp"
#include "blas.hpp"
#ifdef HAVE_CUBLAS
//#warning HAVE_CUBLAS
#include "cublas.hpp"
#include <boost/numeric/bindings/cublas/cublas.hpp>
#endif


int main(int argc, char *argv[]) {

    size_t no = atoi(argv[1]);
    size_t nv =  atoi(argv[2]);
    size_t n = no;

    using namespace cc;
    using tensor::as_matrix;

    size_t nops = nv*(2*no*no*no*nv);
    utility::timer timer;

    Symbol< tensor_reference<4> >::range r(0,n);
    Symbol< tensor_reference<4> >::range ro(0,no);
    Symbol< tensor_reference<4> >::range rv(0,nv);
    Symbol< tensor_reference<4> > S;


    {
	Buffer<double> buffer;

	buffer[0].resize(n*no *nv);
	buffer[1].resize(n*no *nv);
	buffer[2].resize(n*no *nv);

	for (int i = 0; i < 10; ++i) {
	    buffer[3+i].resize(n*no*nv);
	}

	S.set("a", buffer[0], r, ro, rv, 0);
	S.set("b", buffer[1], r, ro, rv, 0);
	S.set("c", buffer[2], r, ro, rv, 0);

// 	timer.reset();
// #pragma omp parallel 
// 	{
// 	    Symbol< tensor_reference<4> > S;
// 	    S.set("a", buffer[0], r, ro, rv, 0);
// 	    S.set("b", buffer[1], r, ro, rv, 0);
// #pragma omp for schedule(dynamic,1)
// 	    for (int i = 0; i < nv; ++i) {
// 		S.set("c", buffer[3+omp::thread()], r, ro, rv, 0);
// 		BOOST_AUTO(c, (as_matrix<1,2>(S["c"][0])));
// 		blas::gemm(1, as_matrix<1,1>(S["a"][0][i]),
// 			   (as_matrix<1,2>(S["b"][0])),
// 			   1, c);
// 	    }
// 	}
// 	std::cout << timer << " rate: " << nops/double(timer) <<std::endl;

	timer.reset();
#pragma omp parallel for schedule(dynamic,1)
	for (int i = 0; i < nv; ++i) {
	    BOOST_AUTO(c, (as_matrix<1,1>(S["c"][0][i])));
	    blas::gemm(1, as_matrix<1,2>(S["a"][0]),
		       trans(as_matrix<1,2>(S["b"][0])),
		       1, c);
	}
	std::cout << timer << " rate: " << nops/double(timer) <<std::endl;


	S.set("a", buffer[0], r, ro, rv, 0);
	S.set("b", buffer[1], r, ro, rv, 0);
	S.set("c", buffer[2], r, rv, ro, 0);

	timer.reset();
#pragma omp parallel for schedule(dynamic,1)
	for (int i = 0; i < nv; ++i) {
	    BOOST_AUTO(c, (as_matrix<2,1>(S["c"][0])));
	    blas::gemm(1, 
		       trans(as_matrix<1,2>(S["a"][0])),
		       (as_matrix<1,1>(S["b"][0][i])),
		       1, c);
	}
	std::cout << timer << " rate: " << nops/double(timer) <<std::endl;

// 	timer.reset();
// #pragma omp parallel for schedule(dynamic,1)
// 	for (int i = 0; i < nv; ++i) {
// 	    BOOST_AUTO(c, (as_matrix<2,1>(S["c"][0])));
// 	    blas::gemm(1, 
// 		       trans(as_matrix<1,1>(S["b"][0][i])),
// 		       (as_matrix<1,2>(S["a"][0])),
// 		       1, c);
// 	}
// 	std::cout << timer << " rate: " << nops/double(timer) <<std::endl;

    }

    {
	cublas::init();

	double *buffer[] = {
	    cublas::alloc<double>(n*no *nv),
	    cublas::alloc<double>(n*no *nv),
	    cublas::alloc<double>(n*no *nv)
	};

	S.set("a", buffer[0], r, ro, rv, 0);
	S.set("b", buffer[1], r, ro, rv, 0);
	S.set("c", buffer[2], r, ro, rv, 0);

	std::vector<cuda::stream> streams(1);
	std::vector<double*> sb(streams.size());
	foreach (cuda::stream &s, streams) {
	    //s.create();
	}

	// timer.reset();
	// //#pragma omp parallel for schedule(dynamic,1)
	// for (int i = 0; i < nv; ++i) {
	//     cublas::set_stream(streams.at(i%streams.size()));
	//     BOOST_AUTO(c, (as_matrix<1,1>(S["c"][0][i])));
	//     blas::gemm(1, as_matrix<1,2>(S["a"][0]),
	// 	       trans(as_matrix<1,2>(S["b"][0])),
	// 	       1, c, blas::DEVICE);
	// }
	// cuda::synchronize();
	// std::cout << timer << " rate: " << nops/double(timer) <<std::endl;

	timer.reset();
	//#pragma omp parallel for schedule(dynamic,1)
	for (int i = 0; i < nv; ++i) {
	  //cublas::set_stream(streams.at(i%streams.size()));
	    //S.set("c", sb.at(i%streams.size()), r, ro, rv, 0);
	    BOOST_AUTO(c, (as_matrix<1,2>(S["c"][0])));
	    streams.at(i%streams.size()).synchronize();
	    blas::gemm(1, as_matrix<1,1>(S["a"][0][i]),
		       (as_matrix<1,2>(S["b"][0])),
		       1, c, blas::DEVICE);
	}
	cuda::synchronize();
	std::cout << timer << " rate: " << nops/double(timer) <<std::endl;


	S.set("a", buffer[0], r, ro, rv, 0);
	S.set("b", buffer[1], r, ro, rv, 0);
	S.set("c", buffer[2], r, rv, ro, 0);

	timer.reset();
	for (int i = 0; i < nv; ++i) {
	    BOOST_AUTO(c, (as_matrix<2,1>(S["c"][0])));
	    blas::gemm(1, 
		       trans(as_matrix<1,2>(S["a"][0])),
		       (as_matrix<1,1>(S["b"][0][i])),
		       1, c,  blas::DEVICE);
	}
	cuda::synchronize();
	std::cout << timer << " rate: " << nops/double(timer) <<std::endl;

	timer.reset();
	for (int i = 0; i < nv; ++i) {
	    BOOST_AUTO(c, (as_matrix<2,1>(S["c"][0])));
	    blas::gemm(1, 
		       trans(as_matrix<1,1>(S["b"][0][i])),
		       (as_matrix<1,2>(S["a"][0])),
		       1, c,  blas::DEVICE);
	}
	cuda::synchronize();
	std::cout << timer << " rate: " << nops/double(timer) <<std::endl;

    }

}

