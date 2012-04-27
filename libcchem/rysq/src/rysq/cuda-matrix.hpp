#ifndef _RYSQ_CUDA_MATRIX_HPP_
#define _RYSQ_CUDA_MATRIX_HPP_

#include "boost/cuda/array.hpp"
#include "rysq/cuda-memory.hpp"
#include "rysq/cuda-memory.hpp"
#include "rysq/fock.hpp"

namespace rysq {
namespace cuda {

	template<typename T, class A_ = boost::cuda::vector<T> >
	struct block_matrix : block_matrix_base<T> {
	    typedef A_ array_type;
	    typedef block_matrix_base<T> base;
	    using base::size1;
	    using base::size2;
	    using base::size;
	    using base::block;
	    block_matrix() {}
	    T* block(int a,int b) {
		return origin() +  base::layout_.block_at(a,b);
	    }
	    const T* block(int a,int b) const {
		return origin() + base::layout_.block_at(a,b);
	    }
	    void clear() {
		fill(data_.begin(), size(), char(0));
	    }
	    template<class M>
	    void assign(const M &A) {
		base::check_size(A);
		rysq::block_matrix<T>
		    B(size1(), size2(), block().size1(), block().size2());
		B.assign(A);
		boost::cuda::copy(B.data(), data_.begin(), B.size());
	    }
	    void resize(size_t size1, size_t size2, size_t block1, size_t block2) {
		base::layout_ = block_matrix_layout(size1, size2, block1, block2);
		base::size1_ = size1;  base::size2_ = size2;
		data_.resize(size());
	    }
	    array_type& data() { return data_; }
	    const array_type& data() const { return data_; }
	    T* origin() { return data_.begin().get(); }
	    const T* origin() const { return data_.begin().get(); }
	private:
	    array_type data_;
	};

}
}


#endif /* __RYSQ_CUDA_MATRIX_HPP_ */
