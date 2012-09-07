#ifndef RYSQ_FOCK_HPP
#define RYSQ_FOCK_HPP

#include "rysq-core.hpp"

namespace rysq {


    struct block_matrix_layout {

	struct block_structure {
	    block_structure() {}
	    block_structure(size_t size1, size_t size2,
			    size_t stride1, size_t stride2)
		:  size1_(size1), size2_(size2),
		   stride1_(stride1), stride2_(stride2) {}
	    const size_t& size1() const { return size1_; }
	    const size_t& size2() const { return size2_; }
	    const size_t& stride1() const { return stride1_; }
	    const size_t& stride2() const { return stride2_; }
	private:
	    size_t size1_, size2_;
	    size_t stride1_, stride2_;
	};

	block_matrix_layout() {}
	block_matrix_layout(size_t size1,size_t size2,
			    size_t block1, size_t block2) {
	    if (size1%block1 || size2%block2) {
		throw std::exception();
	    }
	    size_t stride1 = block1*block2;
	    block_ = block_structure(block1, block2,
				     stride1, stride1*(size1/block1));
	}
    	size_t element_at(int i,int j) const {
	    int a = i/block_.size1(), b = j/block_.size2();
	    return block_at(a,b) + block_element_at(i - a*block_.size1(),
						    j - b*block_.size2());
	}
	size_t block_at(int a,int b) const {
	    return a*block_.stride1() + b*block_.stride2();
	}
	size_t block_element_at(int i, int j) const { return i + j*block_.size1(); }
	const block_structure& block() const { return block_; }
    private:
	block_structure block_;
    };

    template<typename T, class L = block_matrix_layout>
    struct block_matrix_base {
	typedef L layout;
	typedef typename layout::block_structure block_structure;
	block_matrix_base(){}
	block_matrix_base(size_t size1, size_t size2, size_t block1, size_t block2)
	    : size1_(size1), size2_( size2),
	      layout_(size1, size2, block1, block2) {}
	size_t size1() const { return size1_; }
	size_t size2() const { return size2_; }
	size_t size()const { return size1_*size2_; }
	const block_structure& block() const { return layout_.block(); }
    protected:
	template<class M>
	void check_size(const M &m) const {
	    if (size1() != m.size1() || size2() != m.size2()) throw;
	}
	size_t size1_,  size2_;
	layout layout_;
    };

    template<typename T, class A = std::vector<T> >
    struct block_matrix : block_matrix_base<T> {
	typedef block_matrix_base<T> base;
	block_matrix(size_t size1, size_t size2,
		     size_t block1, size_t block2)
	    : base(size1, size2, block1, block2)
	{
	    data_.reserve(this->size());
	}
	template<class M>
	void assign(const M &m) {
	    base:: check_size(m);
	    for (size_t j = 0; j < m.size2(); ++j) {
		for (size_t i = 0; i <  m.size1(); ++i) {
		    (*this)(i,j) = m(i,j);
		}
	    }
	}
	T& operator()(int i,int j) {
	    return data_[base::layout_.element_at(i,j)];
	}
	const T& operator()(int i,int j) const {
	    return data_[base::layout_.element_at(i,j)];
	}
	T* data() { return &data_[0]; }
	const T* data() const { return &data_[0]; }
    protected:
	A data_;
    };



    template<typename T>
    struct block_matrix_adapter : block_matrix_base<T> {
	typedef block_matrix_base<T> base;
	using base::size1;
	using base::size2;
	block_matrix_adapter() : data_(NULL) {}
	block_matrix_adapter(size_t size1, size_t size2, T *data,
			     size_t block1, size_t block2)
	    : base(size1, size2, block1, block2),
	      data_(data) {}
	virtual ~block_matrix_adapter() {}
	template<class M>
	void assign(const M &A) {
	    if (A.size1() != size1() || A.size2() != size2()) throw;
	    for (size_t j = 0; j < size2(); ++j) {
		for (size_t i = 0; i < size1(); ++i) {
		    (*this)(i,j) = A(i,j);
		}
	    }
	}

	T& operator()(int i,int j) {
	    return data_[base::layout_.element_at(i,j)];
	}
	const T& operator()(int i,int j) const {
	    return data_[base::layout_.element_at(i,j)];
	}

	using base::block;

	virtual T* block(int i, int j) {
	    return data_ + base::layout_.block_at(i,j);
	}
	virtual	const T* block(int i, int j) const {
	    return data_ + base::layout_.block_at(i,j);
	}

	void bind(T * origin) { data_ = origin; }
	T* data() { return data_; }
	const T* data() const { return data_; }
    protected:
	T* data_;
	int base_[2];
    };


    namespace hf {

	template<class M>
	struct matrix_set : index_set<M> {
	    matrix_set(const int (&base)[4]) {
		std::copy(base, base + 4, base_index.begin());
	    }
	    boost::array<int,4> base_index;
	};

	template<class M>
	struct matrix_ptr_set : matrix_set<M*> {
	    typedef matrix_set<M*> base_type;
	    matrix_ptr_set(const int (&base)[4]):base_type(base){}
	    template<class M2>
	    matrix_ptr_set(matrix_set<M2> &set)
		: base_type(set.base_index.elems)
	    {
		for (int i = 0; i < 12; ++i)
		    (*this)[i] = &set[i];
	    }
	};

    }

    class Fock {
    public:

	typedef block_matrix_adapter<const double> density_matrix;
	typedef block_matrix_adapter<double> fock_matrix;

	typedef hf::matrix_ptr_set<density_matrix> density_matrix_set;
	typedef hf::matrix_ptr_set<fock_matrix> fock_matrix_set;

	typedef index_set<double*> fock_set;
	typedef index_set<const double*> density_set;

	class Impl;

	struct Parameters {
	    double cscale, xscale, cutoff;
	    Parameters(double cscale, double xscale, double cutoff = 1.0e-10)
		: cscale(cscale), xscale(xscale), cutoff(cutoff) {}
	    void scale(double scale) {
		cscale *= scale;
		xscale *= scale;
	    }
	};

	Fock(const Quartet<Shell> &shells);
	~Fock();

	void operator()(const Quartet<Center> &centers,
			density_set D, fock_set F,
			const Parameters &parameters);

	void operator()(const std::vector<Center> &centers,
			const std::vector<Int4> &quartets,
			density_matrix_set D, fock_matrix_set F,
			const Parameters &parameters);

	static void eval(const Quartet<type> &quartet,
			 const std::vector<Int4> &quartets,
			density_matrix_set D, fock_matrix_set F,
			 const double *Eri, const Parameters &parameters);

	static void eval(const Quartet<type> &quartet,
			 density_set D, fock_set F,
			 const double *Eri, const Parameters &parameters);
    private:
	Transpose transpose_;
	void *kernel_;
    };

}


#endif /* RYSQ_FOCK_HPP */
