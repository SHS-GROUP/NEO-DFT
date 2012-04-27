#ifndef _MATRIX_BASE_HPP_
#define _MATRIX_BASE_HPP_

namespace matrix {

    // using namespace boost;

    typedef size_t size_type;
    typedef int index_type;

    typedef size_type shape_type[2];

    template<typename T>
    class matrix_base {
    public:
	typedef T value_type;
	typedef T element;
	typedef T& reference;
	typedef T* pointer;
	
	const shape_type& shape() const { return shape_; }
	size_type size() const { return size_; }
	size_type size1() const { return shape_[0]; }
	size_type size2() const { return shape_[1]; }


    protected:
	pointer base_;
	shape_type shape_;
	size_type size_;

	void init(const shape_type &shape) {
	    size_ = 1;
	    for (int i = 0; i < 2; ++i) {
		shape_[i] = shape[i];
		size_ *= shape_[i];
	    }
	}

// 	template<int D>
// 	index_range range(const type<dim_index<D,element> >&) const {
// 	    return util::range<index_type>(shape_[D]);
// 	}

// 	template<class C>
// 	index_range range(const type<col_index<C> >&) const {
// 	    typedef dim_index<col_index<C>::dim, C> dim_index;
// 	    return range(type<dim_index>());
// 	}

    };


}

#endif /* _MATRIX_BASE_HPP_ */
