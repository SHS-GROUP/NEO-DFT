#ifndef ADAPTER_RYSQ_HPP
#define ADAPTER_RYSQ_HPP

#include <rysq.hpp>
#include "basis/basis.hpp"
#include "matrix/matrix.hpp"

namespace adapter {
namespace rysq {

    struct Shell : ::rysq::Shell {
	typedef ::rysq::Shell super;
	Shell(const Basis::Shell::Data &shell)
	    : super(::rysq::type(shell.type), shell.K(),
		    shell.a(), &shell.C().front()) {
	}
    };

    struct Quartet {
	Quartet(const Basis::Shell::Data &A, const Basis::Shell::Data &B,
		const Basis::Shell::Data &C, const Basis::Shell::Data &D)
	    : A_(A), B_(B), C_(C), D_(D),
	      quartet_(A_, B_, C_, D_) {}
	typedef ::rysq::Quartet< ::rysq::Shell >& reference;
	operator reference() { return quartet_; }
    private:
	adapter::rysq::Shell A_, B_, C_, D_;
	::rysq::Quartet< ::rysq::Shell > quartet_;
    };

    template<typename T>
    struct block_matrix : ::rysq::block_matrix_adapter<T> {
	typedef ::rysq::block_matrix_adapter<T> base;
	block_matrix() {}
	template<class Matrix>
	block_matrix(Matrix &A)
	    : base(A.size1(), A.size2(), A.data(), A.block1(), A.block2()) {}
    };

    struct fock_matrix : block_matrix<double> {
	typedef block_matrix<double> base;
	fock_matrix() {}
	template<class Matrix>
	fock_matrix(Matrix &A) : base(A) {}
    };

    struct density_matrix : block_matrix<const double> {
	typedef block_matrix<const double> base;
	density_matrix() {}
	template<class Matrix>
	density_matrix(const Matrix &A) : base(A){}
    };

}
}

#endif // ADAPTER_RYSQ_HPP
