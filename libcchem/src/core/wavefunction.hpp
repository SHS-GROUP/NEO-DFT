#ifndef WAVEFUNCTION_HPP
#define WAVEFUNCTION_HPP

#include "basis/basis.hpp"
#include "basis/normalize.hpp"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace detail {
    struct wavefunction {
	typedef Basis basis_type;
    };
}

struct Wavefunction {
    typedef detail::wavefunction::basis_type Basis;
    typedef boost::numeric::ublas::matrix<double,
					  boost::numeric::ublas::column_major>
	matrix_type;
    typedef boost::numeric::ublas::matrix_range<const matrix_type> matrix_range;
    typedef boost::numeric::ublas::vector<double> vector_type;

    struct Orbitals {
	Orbitals(size_t start, size_t stop)
	    : size_(stop - start), start_(start), stop_(stop) {}
	size_t size() const { return size_; }
	size_t start()const {return start_; }
	size_t stop() const { return stop_; }
	operator boost::numeric::ublas::range() const {
	    return boost::numeric::ublas::range(start_, stop_);
	}
    private:
	size_t size_, start_, stop_;
    };

    template<class C, class E>
    Wavefunction(const Basis& basis,
		 size_t core, size_t active, size_t virtuals,
		 const boost::numeric::ublas::matrix_expression<C> &coefficients,
		 const boost::numeric::ublas::vector_expression<E> &energies)
	: basis_(basis),
	  core_(0, core),
	  active_(core_.stop(), core_.stop() + active),
	  virtuals_(active_.stop(), active_.stop() + virtuals),
	  coefficients_(coefficients),
	  energies_(energies) {}

    const Basis& basis() const { return basis_; }
    const Orbitals& core() const { return core_; }
    const Orbitals& active() const { return active_; }
    const Orbitals occupied() const { return Orbitals(0, active_.stop()); }
    const Orbitals& virtuals() const { return virtuals_; }
    const vector_type& energies() const { return energies_; }

    const matrix_type& C() const {return coefficients_; }
    matrix_range C(const Orbitals& r) const {
	namespace ublas = boost::numeric::ublas;
	return ublas::project(C(), ublas::range(0, C().size1()), r);
    }

    void normalize() {
	basis::normalize_matrix1(basis_, coefficients_);
    }

private:
    const Basis& basis_;
    Orbitals core_, active_, virtuals_;
    matrix_type coefficients_;
    vector_type energies_;
};

#endif // WAVEFUNCTION_HPP
