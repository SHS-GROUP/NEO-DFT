#ifndef _MATRIX_BLOCK_MATRIX_HPP_
#define _MATRIX_BLOCK_MATRIX_HPP_

#include <boost/type.hpp>
#include <boost/multi_array.hpp>

#include "matrix/matrix.hpp"
#include "matrix/matrix_base.hpp"

#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <numeric>
#include <cassert> 

// #include "iterator/iterator.hpp"
// #include "util/util.hpp"
// #include "util/array.hpp"

namespace ublas = boost::numeric::ublas;

namespace matrix {

    using namespace boost;

    template<typename T>    
    struct block_matrix_semantics {
	typedef int index_type;
	typedef T* block;
	typedef const T* const_block;
        virtual ~block_matrix_semantics() {}
	virtual block b(index_type i, index_type j) = 0;
	virtual const_block b(index_type i, index_type j) const = 0;
    };

    template<typename T>
    struct block_matrix_base : virtual matrix_semantics<T>,
	virtual block_matrix_semantics<T>
    {
  	typedef matrix_semantics<T> semantics;
	typedef typename semantics::index_type index_type;

	typedef size_t size_type;
	typedef size_type shape_type[2];
	typedef T value_type;
	typedef T element;
	typedef T* pointer;
	typedef const T* const_pointer;
  };

    template<typename T>
    class block_matrix : public virtual block_matrix_base<T> {
    public:
	typedef block_matrix_base<T> base;

	typedef typename base::index_type index_type;
	typedef typename base::size_type size_type;
	typedef typename base::shape_type shape_type;

	typedef typename base::reference reference;
	typedef typename base::const_reference const_reference;

	typedef typename base::pointer pointer;
	typedef typename base::const_pointer const_pointer;

	typedef typename base::block block;
	typedef typename base::const_block const_block;

	typedef boost::array<size_type,2> size_array;

	void resize(const size_type *matrix, const size_type *block) {
	    size_t size = this->size();
	    set(matrix); set_block(block);
	    if (this->size() > size) allocate();
	}

	void resize(size_t size1, size_t size2, size_t block1, size_t block2) {
	    size_t size[2] = { size1, size2 };
	    size_t block[2] = {block1, block2 };
	    resize(size, block);
	}

	void initialize() {
	    shape_type shape = { 0, 0 };
	    shape_type block_shape = { 1, 1 };
	    set(shape);
	    set_block(block_shape);
	    this->raw(NULL);
	}	    

	block_matrix() {
	    initialize();
	}

	block_matrix(size_type m, size_type n) {
	    initialize();
	    size_type shape[] = { m, n };
	    size_type block[] = { 1, 1 };
	    resize(shape, block);
	}

	block_matrix(const shape_type &shape, const shape_type &block_shape) {
	    initialize();
	    size_type s[] = { shape[0], shape[1] };
	    size_type b[] = { block_shape[0], block_shape[1] };
	    resize(s, b);
	}

	block_matrix(const shape_type &shape, const shape_type &block_shape,
		     pointer p) {
	    set(shape);
	    set_block(block_shape);
	    raw(p);
	}

        virtual ~block_matrix() { this->free(); }

	void allocate() {
	    free();
	    data_ = new T[size()];
	    allocated_ = true;
	}

	void free() {
	    if (allocated_) delete[] data_;
	    raw(NULL);
	}

	void raw(pointer raw) {
	    allocated_ = false;
	    data_ = raw;
	}

	pointer raw() { return data_; }

	const_pointer raw() const { return data_; }

	pointer data() { return data_; }

	const_pointer data() const { return data_; }

	 const shape_type& shape() const { return shape_; }

	virtual size_type size1() const { return shape()[0]; }

	virtual size_type size2() const { return shape()[1]; }

	 const size_type size() const { return size1()*size2(); }

	size_t block1()const { return block_shape_[0]; }
	size_t block2()const { return block_shape_[1]; }

	 const shape_type& block_shape() const { return block_shape_; }

	size_type block_size() const { return block_size_; }

	void clear() {
	    std::fill(data(), data() + size(), 0);
	}

 	virtual reference operator()(int i, int j) {
	    return data_[reference_index(i,j)];
	}

 	virtual const_reference operator()(int i, int j) const {
	    return data_[reference_index(i,j)];
	}

	virtual block b(int a, int b) { return block_(a,b); }

	virtual const_block b(int a, int b) const { return block_(a,b); }

	size_t block_ld() const {return ld_block_; }

	void operator+=(const block_matrix &m) {
	    check_range(m);
	    if (block1() == m.block1() && block2() == m.block2()) {
		for (size_t i = 0; i < size(); ++i) {
		    data()[i] += m.data()[i];
		}
	    }
	    else {
		this->operator+=<block_matrix>(m);
	    }
	}

	template<class  M>
	void operator += (const M & m) {
	    check_range(m);
	    for (size_t j = 0; j < size2(); ++j) {
		for (size_t i = 0; i < size1(); ++i) {
		    (*this)(i,j) += m(i,j);
		}
	    }
	}

    protected:
	friend class meta_matrix<block_matrix>;
	friend class block_meta_matrix<block_matrix>;
	friend class boost::multi_array<block_matrix,2>;

	shape_type shape_, block_shape_, block_num_;
	size_type block_size_, ld_block_;

	pointer data_;
	bool allocated_;

	template<class M>
	void check_range(const M & m) const {
	    if ( size1() != m.size1() || size2() != m.size2()) 
		throw std::range_error("invalid dimensions");
	}

	template <typename I>
	static bool check_range(const I &i, const I &b, const I &e) {
	    return (b <= i && i < e);
	}

	void set(const size_type *shape) {
	    for (int i = 0; i < 2; ++i) { shape_[i] = shape[i]; }
	}

	void set_block(const size_type *block_shape) {
	    block_size_ = 1;
	    for (int i = 0; i < 2; ++i) {
		assert(shape_[i]%block_shape[i] == 0);
		block_shape_[i] = block_shape[i];
		block_num_[i] = shape_[i]/block_shape_[i];
		block_size_ *= block_shape_[i];
	    }
	    ld_block_ = block_num_[0]*block_size_;
	}

	block block_(int a, int b) const {
	    assert(check_range<int>(a, 0, block_num_[0]));
	    assert(check_range<int>(b, 0, block_num_[1]));
	    return data_ + block_index(a,b);
	}

	index_type block_index(index_type i, index_type j) const {
	    return i*block_size() + j*ld_block_;
	}

	index_type reference_index(index_type i, index_type j) const {
	    assert(check_range<index_type>(i, 0, shape()[0]));
	    assert(check_range<index_type>(j, 0, shape()[1]));
	    index_type a = i/block_shape_[0];
	    index_type b = j/block_shape_[1];
	    i -= a*block_shape_[0]; 
	    j -= b*block_shape_[1]; 
	    //printf("%i,%i\n", block_index(a,b) + (i + j*block_shape_[0]), size());
	    assert((block_index(a,b) + (i + j*block_shape_[0])) < size());
	    return block_index(a,b) + (i + j*block_shape_[0]);
	}	    


	int block_element(int i, int j, int a, int b) const {
	    i -= a*block_shape()[0]; 
	    j -= b*block_shape()[1]; 
	    return i + j*block_shape()[0];
	}

	

    };

}

#endif /* _MATRIX_BLOCK_MATRIX_HPP_ */
