#ifndef BOOST_PARALLEL_HPP
#define BOOST_PARALLEL_HPP

#include <boost/thread/mutex.hpp>
#include <boost/iterator/iterator_facade.hpp>

namespace boost {

    template<class I, class F>
    struct parallel_iterator : public boost::iterator_facade<
        parallel_iterator<I,F>,
	typename I::value_type,
	boost::forward_traversal_tag,
	typename I::reference>
    {
	parallel_iterator(I begin, I end, F &position)
	    : begin_(begin), end_(end), position_(position)
	{
	    iterator_ = begin_;
	    increment();
	}
	parallel_iterator(I end, F &position)
	    : begin_(end), end_(end), position_(position)
	{
	    iterator_ = begin_;
	}
    private:
	friend class boost::iterator_core_access;
	template<class, class> friend class parallel_iterator;
	void increment() const {
	    size_t next = position_++;
	    while (iterator_ != end_ && iterator_ != (begin_+next)) {
		std::advance(iterator_, 1);
	    }
	}
	template<class J, class Q>
	bool equal(const parallel_iterator<J,Q> &other) const {
	    return (this->iterator_ == other.iterator_);
	}
	typename I::reference dereference() const {
	    BOOST_VERIFY((iterator_ != end_));
	    return *iterator_;
	}
    private:
	I begin_, end_;
	F &position_;
	mutable I iterator_;
    };


    template<class R, class F,
	     class I = parallel_iterator<typename boost::mpl::if_<
					     boost::is_const<R>,
					     typename R::const_iterator,
					     typename R::iterator>::type, F> >
    struct parallel_range {
	typedef I iterator;
	typedef parallel_iterator<typename R::const_iterator, F> const_iterator;
	parallel_range(R &r, F &f) : data_(r), increment_(f) {}
	iterator begin() {
	    return iterator(data_.begin(), data_.end(), increment_);
	}
	iterator end() {
	    return iterator(data_.end(), increment_);
	}
	const_iterator begin() const {
	    return const_iterator(data_.begin(), data_.end(), increment_);
	}
	const_iterator end() const {
	    return const_iterator(data_.end(), increment_);
	}
	R& data() { return data_; }
	const R& data() const { return data_; }
    private:

	template<class T, class U>
	static T min(T a, U b) {
	    return (std::distance<T>(a, b) > 0) ? a : b;
	}

	// void initialize() const {
	//     if (initialized_) {
	// 	initialized_ = true;
	// 	increment_++;
	//     }
	// }
	
	R &data_;
	F &increment_;
    };

}

#endif // BOOST_PARALLEL_HPP
