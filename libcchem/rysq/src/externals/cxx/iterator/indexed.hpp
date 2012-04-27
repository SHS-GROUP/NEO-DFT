#ifndef _ITERATOR_INDEXED_HPP_
#define _ITERATOR_INDEXED_HPP_

#include <boost/multi_array/storage_order.hpp>
#include <boost/range.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/type_traits/add_reference.hpp>
#include <iterator>

#include "cxx/array.hpp"

namespace cxx {
namespace iterator {

    using boost::c_storage_order;
    using boost::fortran_storage_order;
    using boost::array;

    template<size_t N, class Order = boost::c_storage_order>
    struct multi_index {
	typedef boost::array<size_t,N> type;
	multi_index(const size_t (&dims)[N]) {
	    std::copy(dims, dims + N, dims_);
	    index1d_ = 0;
	    update();
	}
	size_t dim(size_t i) const {return dims_[i] ; }
	size_t operator[](size_t i) const { return index_[i]; }
	multi_index& operator++() { ++index1d_; update(); return *this; }
	operator const type&() const { return index_; }
	static type index(size_t index1d, const size_t (&dims)[N],
			  const boost::c_storage_order &order) {
	    type v;
	    size_t r = index1d;
	    for (size_t i = N-1; i > 0; --i) {
		size_t q = r/dims[i];
		v[i] = r - q*dims[i];
		r = q;
	    }
	    v[0] = r;
	    // std::cout << v << index1d  << boost::make_array(dims) << "\n";
	    return v;
	}
	static type index(size_t index1d, const size_t (&dims)[N],
			  const boost::fortran_storage_order &order) {
	    type v;
	    size_t r = index1d;
	    for (size_t i = 0; i < N-1; ++i) {
		size_t q = r/dims[i];
		v[i] = r - q*dims[i];
		r = q;
	    }
	    v[N-1] = r;
	    return v;
	}
    private:
	void update() { index_ = index(index1d_, dims_, Order()); }
	size_t dims_[N];
	type index_;
	size_t index1d_;
    };

    using boost::iterator_range;

    template<size_t N, class F>
    class indexed_iterator
	: public boost::iterator_facade<indexed_iterator<N,F>,
				 typename F::reference,
				 boost::forward_traversal_tag>
    {
    public:
	typedef typename F::storage_order_type storage_order_type;
	typedef array<size_t,N> index_array;
	indexed_iterator(const index_array &dims, F &f, size_t index = 0)
	    : dims_(dims), f_(f), index_(index), store()
	{
	    update();
	}
	const index_array& index() const { return index_n_; }
    private:
	friend class boost::iterator_core_access;
	index_array dims_;	
	F &f_;
	mutable size_t index_;
	mutable index_array index_n_;
	const storage_order_type store;
	void increment() const { ++index_; update(); }
	bool equal(const indexed_iterator &other) const {
	    return index() == other.index();
	}
	typename F::reference dereference() { return f_(index()); }
	typename F::const_reference dereference() const { return f_(index()); }
	void update() const { index_n_ = index_n(index_, dims_, store); }
    };

    template<size_t N, class F>
    iterator_range<indexed_iterator<N,F> >
    indexed_range(const array<size_t,N> &dims, F &f, size_t begin, size_t end,
		  const c_storage_order &order = c_storage_order()) {
	typedef indexed_iterator<N,F> it;
	return make_iterator_range(it(dims, f, begin), it(dims, f, end));
    }


    template<class MA>
    iterator_range<indexed_iterator<MA::dimensionality,MA> >
    indexed_range(MA &A) {
	typedef indexed_iterator<MA::dimensionality,MA> iterator;
	return indexed_range(A.shape(), A, 0, A.num_elements());
    }
    

}
}

#endif /* _ITERATOR_INDEXED_HPP_ */

