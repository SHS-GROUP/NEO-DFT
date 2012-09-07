#ifndef BOOST_PARALLEL_RANGE_HPP
#define BOOST_PARALLEL_RANGE_HPP

#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/thread/mutex.hpp>
#include "externals/boost/parallel/range.hpp"

namespace boost {

    struct parallel {

	template<class C>
	struct atomic {
	    template<class A>
	    atomic(const A &a) : data_(C(a)) {}
	    template<class F>
	    C operator()(F f) {
		boost::mutex::scoped_lock lock(mutex_);
		f(data_);
		return data_;
	    }
	    operator C() {
		boost::mutex::scoped_lock lock(mutex_);
		return data_;
	    }
	protected:
	    boost::mutex mutex_;
	    C data_;
	public:
	    struct increment : atomic<C> {
		explicit increment(const C &initial = C()) : atomic<C>(initial) {}
		C operator++(int) {
		    boost::mutex::scoped_lock lock(this->mutex_);
		    C value = this->data_++;
		    return value;
		}
		C operator++() {
		    boost::mutex::scoped_lock lock(this->mutex_);
		    C value = ++this->data_;
		    return value;
		}
	    };	    
	};


	template<class R>
	parallel_range<R, atomic<size_t>::increment>
	range(R &r) const {
	    return parallel_range<R, atomic<size_t>::increment>(r, increment_);
	}

    private:
	mutable atomic<size_t>::increment increment_;
    };

}

#endif // BOOST_PARALLEL_RANGE_HPP
