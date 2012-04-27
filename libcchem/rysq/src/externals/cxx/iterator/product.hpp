#ifndef _CXX_PRODUCT_HPP_
#define _CXX_PRODUCT_HPP_

#include <iterator>
#include <boost/array.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>
#include "cxx/iterator/indexed.hpp"

namespace cxx {

    namespace iterator {


	template<size_t, class>
	struct product_order;


	template<size_t N>
	struct product_order<N, boost::fortran_storage_order> {
	    static const size_t leading = 0, last = N-1;
	    static const int direction = 1;
	};

	template<size_t N>
	struct product_order<N, boost::c_storage_order> {
	    static const size_t leading = N-1, last = 0;
	    static const int direction = -1;
	};


	template<size_t N, class I, class Enable = void>
	struct product_tuple {
	    typedef typename std::iterator_traits<I>::value_type value_type;
	    typedef typename std::iterator_traits<I>::reference reference;
	    typedef typename std::iterator_traits<I>::pointer pointer;
	    reference operator[](size_t i) const { return *data_[i]; }
	    operator boost::array<value_type,N>() const {
		boost::array<value_type,N> array;
		for (size_t i = 0; i < N; ++i) array[i] = *data_[i];
		return array;
	    }
	protected:
	    template<size_t, class, class> friend class product_iterator;
	    boost::array<pointer,N> data_;
	    void bind(size_t i, reference ref) { data_[i] = &ref; }
	};

	template<size_t N, class I>
	struct product_tuple<N, I,
			     typename boost::disable_if<
				 boost::is_reference<typename I::reference> >::type>
	{
	    typedef typename std::iterator_traits<I>::value_type value_type;
	    typedef typename std::iterator_traits<I>::reference reference;
	    reference operator[](size_t i) const { return data_[i]; }
	    operator boost::array<value_type,N>() const { return data_; }
	protected:
	    template<size_t, class, class> friend class product_iterator;
	    boost::array<reference,N> data_;
	    void bind(size_t i, reference ref) { data_[i] = ref; }
	};

	template<size_t N, class I, class Order = boost::c_storage_order>
	class product_iterator
	    : public boost::iterator_facade<product_iterator<N,I,Order>,
					    product_tuple<N,I>,
					    boost::forward_traversal_tag,
					    product_tuple<N,I> >
	{
	public:
	    product_iterator(const I (&its)[N], const size_t (&dims)[N])
		: index_(dims) {
		for (size_t i = 0; i < N; ++i) { its_[i] = its[i]; }
	    }
	private:
	    friend class boost::iterator_core_access;
	    void increment() const { ++index_; }
	    bool equal(const product_iterator &other) const {
		typedef product_order<N,Order> order;
		const size_t last = order::last;
		bool eq = true;
		for (size_t i = 0; i < N;  ++i) {
		    int j = index_[i] ;
		    j += (index_[i] == 0 && i != last)*index_.dim(i);
		    eq = eq && ((its_[i] + j) == other.its_[i]);
		}
		return eq;
	    }
	    product_tuple<N,I> dereference() const {
		product_tuple<N,I> t;
		for (size_t i = 0; i < N; ++i) {
		    t.bind(i, (its_[i])[index_[i]]);
		}
		return t;
	    }
	    I its_[N];
	    mutable multi_index<N,Order> index_;
	};

	template<size_t N, class I, class Order>
	boost::iterator_range<product_iterator<N,I,Order> >
	product(I (&begin)[N], I (&end)[N], Order order) {
	    typedef product_iterator<N,I,Order> it;
	    size_t dims[N];
	    for (size_t i = 0; i < N; ++i)  { dims[i] = end[i] - begin[i]; }
	    // std::cout << boost:: make_array(dims);
	    return boost::make_iterator_range(it(begin, dims), it(end, dims));
	}

	template<size_t N, class C, class Order>
	boost::iterator_range<
	    product_iterator<N, typename C::const_iterator, Order> >
	product(const C &c, Order order) {
	    typename C::const_iterator begin[N], end[N];
	    for (size_t i = 0; i < N; ++i)  {
	    	begin[i] = c.begin();
	    	end[i] = c.end();
	    }
	    return product(begin, end, order);
	}

	template<size_t N, class C, class Order>
	boost::iterator_range<
	    product_iterator<N, typename C::const_iterator, Order> >
	product(const C (&a)[N], Order order) {
	    typename C::const_iterator begin[N], end[N];
	    for (size_t i = 0; i < N; ++i)  {
		begin[i] = a[i].begin();
		end[i] = a[i].end();
		// std::cout << begin[i] << end[i] << std::endl ;
	    }
	    return product(begin, end, order);
	}

	template<size_t N, class C>
	boost::iterator_range<
	    product_iterator<N, typename C::const_iterator, boost::c_storage_order> >
	product(const C &c) {
	    return product<N>(c, boost::c_storage_order());
	}

	template<size_t N, class C>
	boost::iterator_range<
	    product_iterator<N, typename C::const_iterator, boost::c_storage_order> >
	product(const C (&a)[N]) {
	    return product(a, boost::c_storage_order());
	}

    }
}

#endif /* _CXX_PRODUCT_HPP_ */
