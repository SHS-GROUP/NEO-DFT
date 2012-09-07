#ifndef _BASIS_SHELL_HPP_
#define _BASIS_SHELL_HPP_

#include "basis/config.hpp"

#include <vector>
#include <iostream>
#include <boost/array.hpp>

#include <boost/noncopyable.hpp>

namespace basis {
namespace shell  {

    enum type { BASIS_SHELL_ENUM };

}
}

namespace basis {

    typedef boost::array<double,3> Center;

    /**
       @brief
    */
    struct shell_data {
	friend class Basis;

    public:

	struct key_type {
	    key_type() : data_(0) {}
	    key_type(short value) : data_(value) {}
	    bool operator<(const key_type &key) const {
		return (this->data_ < key.data_);
	    }
	    operator short() const { return data_; } 
	    key_type& operator++() {
		++data_;
		return *this;
	    }
	private:
	    short data_;
	};

	typedef key_type Key;
	Key key() const { return key_; }	

	class Function {
	public:
	    struct Range {
		Range(size_t start, size_t finish)
		    : start_(start), size_(finish-start) {}
		size_t size() const { return size_; }
		template<typename T>
		operator T() const {
		    return T(start_, start_ + size_);
		}
	    private:
		size_t start_, size_;
	    };
	    static double coefficient(int i, int j, int k);
	    int mask;
	    double C;
	    Function(int l, int m, int n) :
		mask(l + (m << 8) + (n << 16)),
		C(coefficient(l, m, n)) {}
	    int operator()(int i) const { return (0xFF & (mask >> (i*8))); }
	};

	typedef std::vector<Function>::const_iterator function_iterator;
	const std::vector<Function>& functions() const { return functions_; }

    public:

	const shell::type type;

	/**
	   @brief Data
	   @param K Number of primitives
	   @param exps Primitive exponents, exps[K]
	   @param nc Number of contractions
	   @param coeff Contraction coefficients, coeff[nc][K]
	   @param L Data angular momentums, L[nc]
	*/
	shell_data(int K, double *exps, double *C, int L);

	/**
	   @brief SP shell constructor
	   @param exps exponents
	   @param Cs S shell coefficients
	   @param Cp P shell coefficients
	*/
	shell_data(int K, double *exps, double *Cs, double *Cp);

	/** @brief Get number of primitives */
	int K() const { return K_; }

	int L() const { return L_; }

	/**
	   @brief maximum angular momentum
	   @return maximum angular momentum
	*/
	int Lmin() const { return Lmin_; }

	/** @brief Get primitive exponents */
	const double* a() const { return &(a_.at(0)); }

	std::vector<const double*> C() const {
	    std::vector<const double*> C;
	    for (int i = 0; i < nc_; ++i) C.push_back(this->C(i));
	    return C;
	}

	/** @brief Get contraction coefficients */
	const double* C(int k) const { return &(C_.at(k).at(0)); }

	/** shell size number of functions in shell */
	size_t size() const { return size_; }

	/**
	   @brief Get a primitive exponent
	   @param exp index of the primitive
	   @return primitive exponent
	*/
	double operator()(int k) const { return a_[k]; }

	int nc() const { return nc_; };

	/**
	   @brief Get contraction coefficient of a shell for a particular primitive
	   @param exp index of the primitive
	   @param cn
	*/
	double operator()(int k, int c) const { return C_[c][k]; }

	bool operator==(const shell_data &shell) const;

    private:
	/** number of primitives */
	int K_;
	/** exponents */
	std::vector<double> a_;
	/** number of contractions */
	int nc_;
	/** contraction coefficients */
	std::vector<std::vector<double> > C_;
	/** angular momentums */
	/** Max shell angular momentum */
	int L_;
	int Lmin_;
	int size_;

	std::vector<Function> functions_;
	Key key_;
    };


    struct Shell {

	typedef shell_data Data;
	typedef Data::Function Function;

	size_t size() const { return this->data_->size(); }

	size_t index() const { return index_; }
	size_t start() const { return start_; }
	size_t stop() const { return start() + size(); }
	Function::Range range() const {return Function::Range(start(), stop()); }

	bool operator>(const Shell &rhs) const {
	    return (this->data_ > rhs.data_);
	}

	bool operator>=(const Shell &rhs) const {
	    return (this->data_->L() >= rhs.data_->L());
	}

	const Data* data() const { return data_; }
	operator const Data&() const { return *data_; }

	const std::vector<Data::Function>& functions() const {
	    return data_->functions();
	}

	const Center& center() const { return center_; }
	const double& center(size_t i) const { return center_[i]; }
	void recenter(const Center &r) { center_ = r; }
	int atom()  const { return atom_; }

	int index(const std::vector<Shell> &shells) const {
	    if (shells.empty() || (&shells.front() > this) || (this > &shells.back())) {
		throw std::range_error("shell out of range");
	    }
	    return (this - &shells.front());
	}			   

    private:

	friend class Basis;
        Shell(const Data *data, const Center &center, int atom)
	    : data_(data), center_(center), atom_(atom),
	      index_(0), start_(0) {}
	
	const Data *data_;
	Center center_;
	int atom_;
	int index_;
	int start_;
    };

    std::ostream& operator<<(std::ostream& os, const shell_data &shell);
    std::ostream& operator<<(std::ostream& os, const Shell &shell);

}


#endif /* _BASIS_SHELL_HPP_ */
