#include "runtime.hpp"

#include "bindings/bindings.h"
#include "bindings/bindings.hpp"

#include "hf/hf.hpp"
#include "basis/normalize.hpp"
#include "matrix/matrix.hpp"
#include "parallel.hpp"

template<class V, class I>
V permute(const V &v, const I &index) {
    V u;
    for (size_t i = 0; i < index.size(); ++i) {
	u.push_back(v[index[i]]);
    }
    return u;
}

int CChem_hf_fock(int basis_handle, const double *D, double *F,
		  double cscale, double xscale,
		  const double *screen, double tolerance) {

    try {

	namespace ublas = boost::numeric::ublas;
	Basis &basis = *bindings::any().find<Basis*>(basis_handle);

	struct {
	    std::vector<int> function, shell;
	} index;
	index.function = basis.index();
	index.shell = basis.shell_index();

	basis.sort();
	index.shell = permute(index.shell, basis.shell_index());
	index.function = permute(index.function, basis.index());

	basis.reverse();
	index.shell = permute(index.shell, basis.shell_index());
	index.function = permute(index.function, basis.index());

	int N = basis.num_shells();
	int n = basis.num_functions();

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
	hf::Matrix Dmax(N, N);
	hf::Matrix Kmax(N, N);

	typedef matrix::permutation<int> P;
	matrix::const_symmetric_adapter<double>  screen_adapter(screen,N);
	matrix::assign(Kmax, P(index.shell)(screen_adapter));

	P Pf(index.function);

	const matrix::const_symmetric_adapter<double> D_adapter(D,n);
	D_ = Pf(D_adapter);
	hf::Dmax(basis, D_, Dmax);
	basis::normalize<2>(basis.N(), D_);

	F_ = 0;

	Parallel pe;
	hf::Parameters parameters(cscale, xscale, tolerance);
	hf::fock(basis, D_, F_, parameters, pe, &Kmax, &Dmax);
	basis::normalize<2>(basis.N(), F_);

	matrix::symmetric_adapter<double> F_adapter(F,n);
	Pf(F_adapter) = F_;

    }
    catch (cchem::exception) { return 0; }

    return 1;
}
