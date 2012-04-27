#ifndef PARALLEL_COUNTER_HPP
#define PARALLEL_COUNTER_HPP

#include <boost/thread/locks.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/noncopyable.hpp>

namespace parallel {

    template<typename T = size_t>
    struct counter : boost::noncopyable {
	counter() : initial_(0) { reset(); }
	// counter(const counter &other) : initial_(other.initial_) { reset(); }
	void reset();
	T next();
	T operator++(int) { return next(); }
    private:
	typedef boost::lock_guard<boost::mutex> lock_guard;
	boost::mutex mutex_;
	T initial_;
   };

}

#endif // PARALLEL_COUNTER_HPP
