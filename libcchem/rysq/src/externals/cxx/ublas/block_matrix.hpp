#ifndef _CXX_UBLAS_BLOCK_MATRIX_H_
#define _CXX_UBLAS_BLOCK_MATRIX_H_

namespace cxx {
    namespace ublas {

	namespace detail {

	    template<typename T>
	    struct scalar {
		typedef T& reference;
		typedef const T& const_reference;
		typedef T* pointer;
		typedef const T* const_pointer;
	    };

	}

	template<typename T, class A>
	struct block_matrix {

	    typedef typename detail::scalar<T>::reference reference;
	    typedef typename detail::scalar<T>::const_reference const_reference;
	    typedef typename detail::scalar<T>::pointer pointer;
	    typedef typename detail::scalar<T>::const_pointer const_pointer;

	    typedef A array_type;

	    block_matrix() : data_() {
		resize(0,0, false);
		block_resize(0,0, false);
	    }

	    reference operator()(int i,int j) {
		return reference::bind(data_ + element_offset(i,j));
	    }
	    const reference operator()(int i,int j) const {
		return const_reference::bind(data_ + element_offset(i,j));
	    }

	    void resize(size_t size1, size_t size2, bool preserve = true) {
		if (preserve) throw "not implemented";
		// assign(size_, size1, size2);
	    }

	    array_type& data() { return data_; }
	    const array_type& data() const { return data_; }

	public:
	    pointer block(int i, int j) {
		return data_ + block_offset(i,j);
	    }
	    const pointer block(int i, int j) const {
		return data_ + block_offset(i,j);
	    }
	    size_t block_size1() { return block_[0]; }
	    size_t block_size2() { return block_[1]; }
	    void block_resize(size_t size1, size_t size2, bool preserve = true) {
		if (preserve) throw "not implemented";
		// assign(size_, size1, size2);
	    }

	protected:
	    size_t block_offset(int a, int b) const {
		return a + b*(block_[0]* block_[1]);
	    }
	    size_t element_offset(int i, int j) const {
		int a = i/block_[0];
		int b = j/block_[1];
		size_t offset = (i - a*block_[0]) + ((j - b*block_[1])*block_[0]);
		return block_offset(a,b) + offset;
	    }
	    array_type data_;
	    size_t size_[2], block_[2];
	};

    }
}

#endif /* _CXX_UBLAS_BLOCK_MATRIX_H_ */

