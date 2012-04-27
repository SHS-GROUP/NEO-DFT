#ifndef _UTIL_REDUCE_HPP_
#define _UTIL_REDUCE_HPP_

#include <algorithm>

namespace util {

    template<typename T, size_t N>
    T multiply(const T (&A)[N]) {
	T v = A[0];
	for (size_t i = 1; i < N; ++i) { v *= A[i]; }
	return v;
    }

    template<typename T>
    struct reduce {

	class base {
	protected:
	    T &value_;
	    base(T &bound) : value_(bound) {}
	    base(T &bound, const T &init) : value_(bound) { value_ = init; }
	public:
	    const T& value() const { return value_; }
	};

	struct max : base {
	    typedef void result_type;
	    max() : base(unbound, T(0)) {}
	    max(T &bound) : base(bound) {}
	    void operator()(const T &o) {
		base::value_ = std::max(base::value_, o);
	    }
	private:
	    T unbound;
	};

    };

}
   
#endif /* _UTIL_REDUCE_HPP_ */
