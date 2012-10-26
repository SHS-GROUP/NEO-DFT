#ifndef CUDA_EXCEPTION_HPP_
#define CUDA_EXCEPTION_HPP_

#include <string>
#include <iostream>
#include <stdexcept>
#include <cuda_runtime.h>

#include <boost/preprocessor/stringize.hpp>

namespace boost {
namespace cuda {

    struct exception : std::runtime_error {
	explicit exception(const char* message)
	    : std::runtime_error(message) {}
 // {
 // 	    std::cout << message << std::endl;
 // 	    what_ = message;
 // 	}
 // 	virtual const char *what() const throw() { return what_; }
 // 	operator std::exception() const {
 // 	    return std::runtime_error(what_);
 // 	}
 //    private:
 // 	const char* what_;
    };

    inline std::ostream& operator<<(std::ostream& o, const exception &e) {
	return o << e.what();
    }

    struct configuration_error : exception {
	configuration_error(const char* message)
	    : exception(message) {}
    };

   

    inline void throw_(cudaError_t status,
		       const char *file = 0, const char *line = 0) {
    	if (status == cudaSuccess) return;

	std::string what;
	if (file) {
	    what += std::string(file);
	    if (line) what += std::string(":") + line;
	    what += ": ";
	}
	what += cudaGetErrorString(status);

	if (status == cudaErrorInvalidConfiguration)
	    throw configuration_error(what.c_str());
	else
	    throw exception(what.c_str());
    }

    static void check_status() {
	throw_(cudaGetLastError());
    }

#define BOOST_CUDA_CHECK(F)							\
    if ((F) != cudaSuccess) boost::cuda::throw_					\
				( cudaGetLastError(),				\
				  __FILE__, BOOST_PP_STRINGIZE(__LINE__ ) )

}
}

#endif /* _CUDA_EXCEPTION_HPP_ */
