#ifndef RYSQ_ASSERT_HPP
#define RYSQ_ASSERT_HPP

#include <string>
#include <stdexcept>
#include <sstream>
#include <boost/current_function.hpp>

namespace rysq {

    struct assertion_failed : std::runtime_error {
	static void throw_(char const *file, long line,
			   char const *function,
			   char const *expr) {
	    std::stringstream ss;
	    ss << file << ":" << line;
	    ss << ": " << function;
	    ss << ": " << expr;
	    throw assertion_failed(ss.str());
	}
    private:
	assertion_failed(const std::string &what)
	    : std::runtime_error(what) {}
    };

}

#define RYSQ_ASSERT(expr)						\
    ((expr) ? ((void)0) : ::rysq::assertion_failed::throw_		\
     (__FILE__, __LINE__, BOOST_CURRENT_FUNCTION, #expr))

#endif // RYSQ_ASSERT_HPP
