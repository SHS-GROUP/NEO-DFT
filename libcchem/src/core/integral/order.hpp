#ifndef INTEGRAL_ORDER_HPP
#define INTEGRAL_ORDER_HPP

#include "core/integral.hpp"

#include <boost/ref.hpp>
#include <boost/array.hpp>
#include <boost/assert.hpp>

struct Integral::Order {
    Order() {
	initialize(0,1,2,3);
    }
    Order(size_t i, size_t j, size_t k, size_t l) {
	initialize(i, j, k, l);
    }

    template<typename T>
    boost::array<boost::reference_wrapper<const T>,4>
    braket(const T &A, const T &B, const T &C, const T &D) const {
	return braket(A,B,C,D);
    }

    template<typename T>
    boost::array<boost::reference_wrapper<T>,4>
    braket(T &A, T &B, T &C, T &D) {
	typedef boost::array<boost::reference_wrapper<T>,4> array_type;
	using boost::ref;
	array_type Q = {{ ref(A), ref(B), ref(C), ref(D) }};
	array_type array = {{ Q[data_[0]], Q[data_[1]], Q[data_[2]], Q[data_[3]] }};
	return array;
    }

    bool switch12() const { return  switch12_; }
    bool switch13() const { return  switch13_; }

    
    void change(T &A, T &B, T &C) const {
	if (switch12_) std::swap(A,B);
	if (switch13_) std::swap(A,C);
    }
    void change(T &A, T &B, T &C) const {
	if (switch12_) std::swap(A,B);
	if (switch13_) std::swap(A,C);
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


#endif // INTEGRAL_ORDER_HPP
