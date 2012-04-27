#ifndef _CUDA_GRID_HPP_
#define _CUDA_GRID_HPP_

#include <vector>
#include <algorithm>

namespace cuda {

    template<typename T, size_t N>
    struct grid {
	typedef T* pointer;
	typedef const T* const_pointer;
	typedef size_t index_type[N];
	typedef size_t size_type[N];
	typedef std::vector<size_t> index_vector;

	struct cell {
	    cell(pointer origin, index_type base, size_type shape) {
		size_ = 1;
		for (int i = 0; i < N; ++i) {
		    base_[i] = base[i];
		    shape_[i] = shape[i];
		    size_ *= shape_[i];
		}
	    }
	    size_t size() const { return size_; }
	    const index_type& base() const { return base_; }
	    const size_type& shape() const { return shape_; }
	    operator pointer() { return origin_; }    
	    operator const_pointer() const { return origin_; }
	private:
	    pointer origin_;
	    index_type base_;
	    size_type shape_;
	    size_t size_;
	};
 
	grid(const std::vector<size_t> &dimensions) {}
	size_t size() const { return size_; }
	void bind(pointer origin) { origin_ = origin; }
	cell find(const index_type &index) {
	    index_type base;
	    size_type size;
	    lower_bound(index, base, size);
	    return cell(origin_ + project(base), base, size);
	}
	void lower_bound(const index_type &index, index_type &base,
			 size_type &size) const {
	    for (int i = 0; i < N; ++i) {
		if (index[i] >= index_[i].back()) throw;
		index_vector::const_iterator it =
		    std::upper_bound(index_[i].begin(), index_[i].end(), index[i]);
		base[i] = *it;
		size[i] = *(it + 1) - *it;
	    }
	}
	size_t project(const index_type &index) const {
	    size_t offset = 0, m = 1;
	    for (int i = 0; i < N; ++i) {
		offset += index[i]*m;
		m *= index_[i].back();
	    }
	    return offset;
	}
    private:
	pointer origin_;
	index_vector index_[N];
	size_t size_;
    };

}

#endif /* _CUDA_GRID_HPP_ */
