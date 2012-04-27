#ifndef _UTIL_GENERATOR_HPP_
#define _UTIL_GENERATOR_HPP_

namespace util {

    template<typename T>
    struct generator {
	T value_;
	const T increment_;
	generator(T first, T increment = 1) : value_(first), increment_(increment) {}
	T operator()() {
	    value_ += increment_;
	    return value_ - increment_;
	}
    };

}

#endif /* _UTIL_GENERATOR_HPP_ */
