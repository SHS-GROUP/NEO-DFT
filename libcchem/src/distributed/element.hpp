#ifndef _DISTRIBUTED_ELEMENT_HPP_
#define _DISTRIBUTED_ELEMENT_HPP_

namespace distributed {

    template<typename T, size_t N>
    struct element {
	typedef distributed::array<T,N> DArray;
	typedef typename DArray::index_array index_array;
	DArray &D;
	index_array index, ld;
	element(DArray &D, const index_array &index) : D(D), index(index), ld(1) {}
	operator T() const {
	    T value;
	    D.get(index, index, &value, ld);
	    return value;
	}
	T operator=(T value) {
	    D.put(index, index, &value, ld);
	    return value;
	}
    };

    template<typename T, size_t N>
    struct const_element {
	typedef distributed::array<T,N> DArray;
	typedef typename DArray::index_array index_array;
	const DArray &D;
	mutable index_array index, ld;
	const_element(const DArray &D, const index_array &index)
	    : D(D), index(index), ld(1) {}
	operator T() const {
	    T value;
	    D.get(index, index, &value, ld);
	    return value;
	}
    };

}

#endif /* _DISTRIBUTED_ELEMENT_HPP_ */
