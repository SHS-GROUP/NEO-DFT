#ifndef RYSQ_CUDA_KERNEL_REDUCE_HPP
#define RYSQ_CUDA_KERNEL_REDUCE_HPP

#include "rysq/cuda.hpp"
#include "cuda/matrix.hpp"
#include "cuda/detail.hpp"

namespace rysq {
    namespace cuda {

	// void reduce(const rysq::matrix_data_array<double> A,
	// 	    double scale, rysq::matrix<double> B);

	// static void reduce(const rysq::matrix_data_array<double> A,
	// 	    double scale, rysq::matrix<double> B) {
	//     reduce <double>(A, scale, B);
	// }

	// static void reduce(const rysq::matrix_data_array<double> A,
	// 		   double scale, rysq::matrix<double> B) {
	//     reduce(A, scale, B);
	// }

	void reduce(const cuda::Quartets &quartets,
		    const cuda::detail::Fock::fock_set &fock,
		    cuda::Fock::fock_matrix_set &F,
		    const double (&scale)[6]);

	void reduce(const rysq::matrix_data_array<double> A,
		    double scale, rysq::matrix<double> B,
		    cudaStream_t stream = 0);

    }
}

#endif // RYSQ_CUDA_KERNEL_REDUCE_HPP
