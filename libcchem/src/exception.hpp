#ifndef EXCEPTION_HPP
#define EXCEPTION_HPP

#include <stdio.h>
#include <stdarg.h>

#include <string>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <boost/current_function.hpp>

namespace cchem {

    struct exception : std::runtime_error {
	static std::string str(char const *file, long line,
			       char const *arg0 = NULL,
			       char const *arg1 = NULL) {
	    std::stringstream ss;
	    ss << file << ":" << line;
	    if (arg0) ss << ": " << arg0;
	    if (arg1) ss << ": " << arg1;
	    return ss.str();
	}
	exception(const std::string &what = "")
	    : std::runtime_error(what) {}
	exception(const char* file, int line)
	    : std::runtime_error(str(file, line)) {}
	exception(const char* file, int line, const char* what)
	    : std::runtime_error(str(file, line, what)) {}
	exception(char const *file, long line,
		  char const *function,
		  const std::string &what)
	    : std::runtime_error(str(file, line, function, what.c_str())) {}
    };

    struct assertion_failed : std::runtime_error {
	static void throw_(char const *file, long line,
			   char const *function,
			   char const *expr) {
	    throw assertion_failed
		(cchem::exception::str(file, line, function, expr));
	}
    private:
	assertion_failed(const std::string &what)
	    : std::runtime_error(what) {}
    };

    inline void message(char const *file, long line, const char *format, ...) {
	FILE * stream = stdout;
	fprintf(stream, "%s:%ld: ", file, line);
	va_list args;
	va_start(args, format);
	vfprintf(stream, format, args);
	va_end(args);
    }

    inline void printf(FILE *stream, const char *format, ...) {
	va_list args;
	va_start(args, format);
	vfprintf(stream, format, args);
	va_end(args);
    }

    inline void error(char const *file, long line, const char *format, ...) {
	FILE * stream = stderr;
	fprintf(stream, "%s:%ld: ", file, line);
	va_list args;
	va_start(args, format);
	vfprintf(stream, format, args);
	va_end(args);
	throw exception("");
    }

}

#define CCHEM_EXCEPTION(what)					\
    cchem::exception						\
    (__FILE__, __LINE__, BOOST_CURRENT_FUNCTION, (what))

#define CCHEM_ASSERT(expr)						\
    ((expr) ? ((void)0) : ::cchem::assertion_failed::throw_		\
     (__FILE__, __LINE__, BOOST_CURRENT_FUNCTION, #expr))

#define CCHEM_ERROR(...)			\
    ::cchem::error(__FILE__, __LINE__, __VA_ARGS__)

#define CCHEM_STDOUT(...)			\
    ::cchem::printf(stdout, __VA_ARGS__)

#define CCHEM_MESSAGE(...)			\
    ::cchem::message(__FILE__, __LINE__, __VA_ARGS__)

#endif // EXCEPTION_HPP
