#ifndef UTILITY_ITERATOR_RANGE_HPP
#define UTILITY_ITERATOR_RANGE_HPP


#include <boost/range/iterator_range.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <cassert>
#include "utility/iterator/increment.hpp"

namespace utility {
namespace iterator {

    template <typename T>
    boost::iterator_range<increment_iterator<T> > range(T from, T to, T increment = 1) {
	assert((increment >= T() && from <= to) || (increment < T() && from >= to));
	typedef increment_iterator<T> iterator;
	return boost::make_iterator_range(iterator(from, increment), iterator(to));
    }

    template <typename T>
    boost::iterator_range<increment_iterator<T> > range(T to) {
	typedef increment_iterator<T> iterator;
	return boost::make_iterator_range(iterator(T(0)), iterator(to));
    }

    template<typename T = int>
    struct Range {
	typedef increment_iterator<T> iterator;
	typedef const increment_iterator<T> const_iterator;
	int first_, last_, increment_;
	Range(int last) : first_(0), last_(last), increment_(1) {}
	Range(int first, int last, int increment = 1)
	    : first_(first), last_(last), increment_(increment) {}
	iterator begin() const { return iterator(first_, increment_); }
	iterator end() const { return iterator(last_, increment_); }
    };

}
}

#endif /* UTILITY_ITERATOR_RANGE_HPP */
