#ifndef _ITERATOR_INCREMENT_HPP_
#define _ITERATOR_INCREMENT_HPP_

#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <cassert>

namespace iterator {

    namespace details {

	template<typename T>
	struct increment_cast { typedef T type; };

	template<typename T>
	struct increment_cast<T*> { typedef size_t type; };

    }

    template <typename T>
    class increment_iterator
	: public boost::iterator_facade<increment_iterator<T>,
					const T,
					std::forward_iterator_tag>
    {
    private:
	friend class boost::iterator_core_access;
	typedef typename details::increment_cast<T>::type cast;
	T value, incr;
	void increment() { value += cast(incr); }
	bool equal(const increment_iterator& other) const {
	    //this is probably somewhat problematic, assuming that the "end iterator"
	    //is always the right-hand value?
	    return ((incr >= 0 && value >= other.value) ||
		    (incr < 0 && value <= other.value));
	}
	const T& dereference() const { return value; }
    public:
	increment_iterator(T value, T incr = T(1)): value(value), incr(incr) {}
    };

}


#endif /* _ITERATOR_INCREMENT_HPP_ */
