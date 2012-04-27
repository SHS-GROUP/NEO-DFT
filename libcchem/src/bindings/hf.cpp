#include "bindings/bindings.h"
#include "bindings/bindings.hpp"

#include "hf/hf.hpp"
#include "basis/normalize.hpp"
#include "util/util.hpp"
#include "matrix/matrix.hpp"

#include "parallel/parallel.hpp"

void HF_fock(const int basis_handle, const double *D, double *F,
	     double cscale, double xscale,
	     const double *screen, double tolerance) {
    parallel::counter<size_t> counter;
    if (parallel::Context().smp().rank() != 0) return;

    namespace ublas = boost::numeric::ublas;
    const Basis &basis = *bindings::any().find<Basis*>(basis_handle);
    int N = basis.num_shells();
    //int N2 = binomial2 (N);
    int n = basis.num_functions();
    // int n2 = binomial2 (n);

    typedef matrix::symmetric_adapter<double> symmetric;
    typedef matrix::const_symmetric_adapter<double> const_symmetric;
    typedef matrix::block_meta_matrix<BlockMatrix> MetaMatrix;
    typedef MetaMatrix::size_vector size_vector;

    size_vector block_sizes, shell_sizes, matrix_sizes;
    foreach (const Basis::Block& block, basis.blocks()) {
	shell_sizes.push_back(block.shell().size());
	block_sizes.push_back(block.size());
	// std::cout << block.shell().size() << " "
	// 	      << block.size() << std::endl;
	matrix_sizes.push_back(shell_sizes.back()*block_sizes.back());
    }
    const size_vector matrix_dims[] = { matrix_sizes, matrix_sizes };
    const size_vector matrix_block_dims[] = { shell_sizes, shell_sizes };
    MetaMatrix D_(matrix_dims, matrix_block_dims);
    MetaMatrix F_(matrix_dims, matrix_block_dims);

    // matrix::meta_matrix<Matrix> Dmax(block_sizes);//, Kmax(block_sizes);
    Matrix Dmax(N, N);

    // Matrix K_(N,N);
    // typedef matrix::meta_matrix<Matrix> meta_matrix;
    // // matrix::meta_matrix_decorator<Matrix, Matrix> Kmax(block_sizes, K_);
    // matrix::meta_matrix<Matrix> Kmax(block_sizes);
    // Kmax = K_;
    Matrix Kmax(N, N);

    typedef matrix::permutation<int> P;
    matrix::const_symmetric_adapter<double>  screen_adapter(screen,N);
    matrix::assign(Kmax, P(basis.shell_permutations())(screen_adapter));

    P Pf(basis.function_permutations());

    const matrix::const_symmetric_adapter<double> D_adapter(D,n);
    D_ = Pf(D_adapter);
    hf::Dmax(basis, D_, Dmax);
    basis::normalize_matrix(basis, D_);

    //matrix::zero(F_);
    F_ = 0;

    //matrix::meta_matrix_decorator<BlockMatrix, MetaMatrix&> d(dims, dims, D_, 6);

    hf::Parameters parameters(cscale, xscale, tolerance);
    hf::fock(basis, D_, F_, parameters, counter, &Kmax, &Dmax);
    basis::normalize_matrix(basis, F_);

    matrix::symmetric_adapter<double> F_adapter(F,n);
    Pf(F_adapter) = F_;

}
