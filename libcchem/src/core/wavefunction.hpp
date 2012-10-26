#ifndef WAVEFUNCTION_HPP
#define WAVEFUNCTION_HPP

#include "basis/basis.hpp"
#include "basis/normalize.hpp"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace detail {
    struct wavefunction {
	typedef Basis basis_type;
    };
}

struct Wavefunction {
    typedef detail::wavefunction::basis_type Basis;
    typedef boost::numeric::ublas::matrix<
	double, boost::numeric::ublas::column_major> matrix_type;
    typedef boost::numeric::ublas::matrix_range<const matrix_type> matrix_range;
    typedef boost::numeric::ublas::vector<double> vector_type;
    typedef boost::numeric::ublas::vector_range<const vector_type> vector_range;

    struct Orbitals {
	Orbitals() : size_(0), start_(0), stop_(0) {}
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

    Wavefunction(const Basis& basis) : basis_(basis) {}

    void set_occupied(size_t start, size_t stop, size_t core) {
	occupied_ = Orbitals(start, stop);
	active_ = Orbitals(start+core, stop);
    }

    void set_virtuals(size_t start, size_t stop) {
	virtuals_ = Orbitals(start, stop);
    }

    template<class E>
    void set_C(const boost::numeric::ublas::matrix_expression<E> &C) {
	C_ = C;
    }

    template<class E>
    void set_e(const boost::numeric::ublas::vector_expression<E> &e) {
	e_ = e;
    }

    template<class E>
    void set_F(const boost::numeric::ublas::matrix_expression<E> &F) {
	F_ = F;
    }

    const Basis& basis() const { return basis_; }

    const Orbitals& occupied() const { return occupied_; }
    const Orbitals& active() const { return active_; }
    const Orbitals& virtuals() const { return virtuals_; }
    Orbitals orbitals() const {
	int start = std::min(occupied_.start(), virtuals_.start());
	int finish = std::max(occupied_.stop(), virtuals_.stop());
	return Orbitals(start, finish);
    }

    const vector_type& e() const { return e_; }
    vector_range e(const Orbitals& r) const {
	return boost::numeric::ublas::project(e_, r);
    }

    const matrix_type& F() const { return F_; }
    const matrix_type& C() const { return C_; }

    matrix_range C(const Orbitals& r) const {
	namespace ublas = boost::numeric::ublas;
	return ublas::project(C(), ublas::range(0, C().size1()), r);
    }

    void normalize() {
	basis::normalize<0>(basis_.N(), C_);
    }

    void sort() {
	basis_.sort();
	reorder();
    }

    void reverse() {
	basis_.reverse();
	reorder();
    }

private:
    Basis basis_;
    Orbitals occupied_, active_, virtuals_;
    matrix_type C_, F_;
    vector_type e_;

    void reorder() {
	vector_type c(basis_.size());
	for (size_t j = 0; j < C_.size2(); ++j) {
	    for (size_t i = 0; i < basis_.size(); ++i) {
		c[i] = C_(basis_.index()[i], j);
	    }
	    column(C_,j).assign(c);
	}
    }

};

#endif // WAVEFUNCTION_HPP
