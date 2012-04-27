#ifndef TENSOR_INDEX_HPP
#define TENSOR_INDEX_HPP

#include <boost/array.hpp>
#include "boost/preprocessor/char.hpp"
#include <boost/preprocessor/seq/for_each.hpp>

#include <boost/mpl/bool.hpp>
#include <boost/mpl/count_if.hpp>
#include <boost/mpl/transform.hpp>

#include <boost/fusion/mpl.hpp>
#include <boost/fusion/container.hpp>
#include <boost/fusion/algorithm.hpp>

#include "boost/fusion/merge.hpp"
#include "boost/fusion/intersection.hpp"
#include "boost/fusion/keys.hpp"

#include <boost/range/iterator_range.hpp>
#include <boost/iterator/counting_iterator.hpp>

namespace tensor {
namespace detail {

    namespace fusion = boost::fusion;
    namespace result_of = boost::fusion::result_of;

    using boost::fusion::map;
    using boost::fusion::merge;
    using boost::fusion::intersection;
    using boost::fusion::keys;

    struct index_range
	: boost::iterator_range<boost::counting_iterator<int> >
    {
	typedef boost::iterator_range<boost::counting_iterator<int> > base;
	typedef base::value_type value_type;
	typedef base::const_iterator const_iterator;
	typedef base::difference_type difference_type;

	struct increment_type {
	    operator difference_type() const { return 1; }
	};
	
	static const_iterator min() {
	    return std::numeric_limits<value_type>::min();
	}

	static const_iterator max() {
	    return std::numeric_limits<value_type>::max();
	}

	index_range() : base(min(), max()) {}
	explicit index_range(const base &range) : base(range) {}
	explicit index_range(int size) : base(0, size) {}
	index_range(int begin, int end) : base(begin, end) {}
	index_range(const_iterator begin, const_iterator end,
		    increment_type increment)
	    : base(begin, end) {}
	
	increment_type increment() const { return increment_type(); }

	index_range operator()(int size) const;
	index_range operator()(int begin, int end) const;

	bool operator!=(const index_range &r) const {
	    return !(static_cast<const base&>(*this) ==
		     static_cast<const base&>(r));
	}
    };


}
}

namespace tensor {


    template<int I>
    struct index_key : boost::mpl::int_<I> {};

    template<int I, class R = detail::index_range>
    struct index : boost::fusion::pair<index_key<I>, R>
    {
	static const int value = I;
	typedef index_key<I> key_type;
	typedef R range_type;
	typedef boost::fusion::pair<key_type, range_type> base;
	typedef typename range_type::value_type value_type;
	typedef typename range_type::const_iterator const_iterator;
	typedef typename range_type::increment_type increment_type;

	// typedef boost::fusion::pair<key_type, range_type> pair;
		       
	// typedef key_type first_type;
	// typedef R second_type;
	
	static const_iterator min() { return range_type::min(); }
	static const_iterator max() { return range_type::max(); }

	index() : base(range_type()), all_(true) {}
	explicit index(const range_type &range) : base(range) {}
	explicit index(int size) : base(range_type(size)) {}
	index(int begin, int end) : base(range_type(begin, end)) {}
	index(const_iterator begin, const_iterator end,
	      increment_type increment)
	    : base(range_type(begin, end, increment)) {}
	index operator=(const index &r);

	const range_type& data() const { return base::second; }

	const_iterator begin() const {
	    return data().begin();
	}
	const_iterator end() const {
	    return data().end();
	}
	increment_type increment() const {
	    return data().increment();
	}

	range_type operator()(int size) const;
	range_type operator()(int begin, int end) const;

	index zero_based() const;
	operator const range_type&() const { return data(); }
	bool operator!=(const index &r) const {
	    return !((r.all_ == all_) && (r.data() != data()));
	}
	bool all() const { return all_; }
    private:
	struct bool_ {
	    explicit bool_(bool value = false) : value(value) {}
	    operator const bool&() const { return value; } 
	    bool value;
	} all_;
    };


    namespace detail {

	template<class T>
	struct is_index
	    : boost::mpl::false_ {};

	template<int I, class R>
	struct is_index<const index<I,R> >
	    : boost::mpl::true_ {};

	template<int I, class R>
	struct is_index<index<I,R> >
	    : boost::mpl::true_ {};


    }


    template<class I>
    struct indices : I {
	template<class T>
	struct is_index :
	    detail::is_index<typename boost::remove_reference<T>::type> {}; 
	static const size_t rank =
	    boost::mpl::count_if<I, is_index<boost::mpl::_> >::value;
	indices(const I &i) : I(i) {}
    };


    namespace index_names {

#define TENSOR_INDEX_NAME(Z,TEXT,V)			\
	static const index<BOOST_PP_CHAR(V)> V =	\
	    index<BOOST_PP_CHAR(V)>();

#define TENSOR_INDEX_NAMES				\
	(a)(b)(c)(d)(e)(f)(g)(h)(i)(j)(k)(l)(m)(n)
 
	BOOST_PP_SEQ_FOR_EACH(TENSOR_INDEX_NAME, nil, TENSOR_INDEX_NAMES);


    }

}

#endif // TENSOR_INDEX_HPP
