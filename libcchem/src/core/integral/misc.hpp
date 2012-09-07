#ifndef INTEGRAL_AUX_HPP
#define INTEGRAL_AUX_HPP

#include "basis/basis.hpp"
#include "foreach.hpp"
#include <boost/multi_array.hpp>
#include <boost/ref.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/assert.hpp>

namespace integral {


struct order {
    order() {
	initialize(0,1,2,3);
    }
    order(size_t i, size_t j, size_t k, size_t l) {
	initialize(i, j, k, l);
    }
    bool switch12() const { return  switch12_; }
    bool switch13() const { return  switch13_; }

    template<typename T>    
    void change(T &A, T &B, T &C) const {
	if (switch12_) std::swap(A,B);
	if (switch13_) std::swap(A,C);
    }
    template<typename T>
    void change(boost::array<T,3> &A) const {
	change(A[0], A[1], A[2]);
    }

private:
    boost::array<size_t,4> data_, braket_;
    bool switch12_, switch13_;
    void initialize(size_t i, size_t j, size_t k, size_t l) {
	boost::array<size_t,4> data = {{ i, j, k, l }};

	// check index range
	size_t sum = 0;
	foreach (size_t o, data) {
	    if (o > 4) throw std::range_error("");
	    sum += o;
	}
	if (sum != 6) throw std::range_error("");

	data_ = data;

	switch12_ = false;
	switch12_ |= (i == 0 && k == 1);
	switch12_ |= (i == 1 && k == 0);
	switch12_ |= (i == 2 && k == 3);
	switch12_ |= (i == 3 && k == 2);

	switch13_ = false;
	switch13_ |= (i == 0 && l == 1);
	switch13_ |= (i == 1 && l == 0);
	switch13_ |= (i == 2 && l == 3);
	switch13_ |= (i == 3 && l == 2);

	BOOST_VERIFY(!(switch12_ && switch13_));

	braket_ = data;
	if (switch12_) std::swap(braket_[1], braket_[2]);
	if (switch13_) std::swap(braket_[1], braket_[3]);
    }
};

struct array : boost::multi_array_ref<double,4> {
    typedef boost::multi_array_ref<double,4> base_type;
    typedef base_type::element value_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    array() : base_type(NULL, boost::extents[0][0][0][0])
    {
	// std::cout << "create" << std::endl;
    }
    array(size_t n1, size_t n2, size_t n3, size_t n4)
	: base_type(NULL, boost::extents[0][0][0][0])
    {
	//std::cout << "create" << std::endl;
	resize(n1, n2, n3, n4);
    }
    void resize(size_t n1, size_t n2, size_t n3, size_t n4) {
	boost::array<size_t,4> shape = {{ n4, n3, n2, n1 }};
	size_t size = n1*n2*n3*n4;
	data_.clear();
	data_.resize(size, 0);
	// update base_type parent
	set_base_ptr(&data_[0]);
	num_elements_ = size;
	reshape(shape);
    }
    size_t size() const { return data_.size(); }
    size_t size1() const { return shape()[3]; }
    size_t size2() const { return shape()[2]; }
    size_t size3() const { return shape()[1]; }
    size_t size4() const { return shape()[0]; }
    double& operator()(int i, int j, int k, int l) {
	return (*this)[l][k][j][i];
    }
    const double& operator()(int i, int j, int k, int l) const {
	return (*this)[l][k][j][i];
    }
    array& operator=(double value) {
	std::fill(data_.begin(), data_.end(), value);
	return *this;
    }
private:
    std::vector<double> data_;
};


struct screening {
    typedef boost::numeric::ublas::matrix<double> matrix;
    screening() : cutoff_(0) {}
    screening(const matrix &K, double cutoff = 1e-10)
	: K_(K), cutoff_(cutoff) {}
    bool test(double value) const {return (std::abs(value) > cutoff_); }
    bool test(int i, int j, int k, int l) const {
	if (K_.size1() <= 0) return true;
	return (K_(i,j)*K_(k,l) > cutoff_);
    }
private:
    const matrix K_;
    double cutoff_;
};


}

#endif // INTEGRAL_AUX_HPP
