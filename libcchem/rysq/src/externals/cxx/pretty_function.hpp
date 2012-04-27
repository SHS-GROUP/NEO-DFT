#ifndef _CXX_PRETTY_FUNCTION_HPP_
#define _CXX_PRETTY_FUNCTION_HPP_

#include <string>
#include <sstream>

namespace pretty_function_ {

    struct string;

    template<typename T>
    string operator+(const string &s, const T &t);

    struct string {
	std::string data_;
	explicit string(const char *function) : data_(function) {}
	explicit string(const std::string &function) : data_(function) {}
	operator std::string&() { return data_; }
	operator const std::string&() const { return data_; }
	string operator()() const { return *this; }
	template<typename T0>
	string operator()(const T0 &t0) const {
	    return (*this + ": " + t0);
	}
	template<typename T0, typename T1>
	string operator()(const T0 &t0, const T1 &t1) const {
	    return (*this)(t0) + t1;
	}
	template<typename T0, typename T1, typename T2>
	string operator()(const T0 &t0, const T1 &t1, const T2 &t2) const {
	    return (*this)(t0) + t1 + t2;
	}
    };

    template<typename T>
    string operator+(const string &s, const T &t) {
    	std::stringstream ss;
    	ss << std::string(s) << t;
    	return string(ss.str());
    }

}

#define PRETTY_FUNCTION(...)						\
    (pretty_function_::string(__PRETTY_FUNCTION__)(__VA_ARGS__))

#endif /* _CXX_PRETTY_FUNCTION_HPP_ */
