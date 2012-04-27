#ifndef boost_cublasbindings_cublasbindings_hpp
#define boost_cublasbindings_cublasbindings_hpp

#include <boost/numeric/bindings/cublas/forward.hpp>
#include <boost/numeric/bindings/cublas/vector.hpp>
#include <boost/numeric/bindings/cublas/matrix.hpp>
#include <boost/numeric/bindings/cublas/host.hpp>
#include <boost/numeric/bindings/cublas/expression.hpp>
#include <boost/numeric/bindings/cublas/exception.hpp>

#include <boost/numeric/bindings/cublas/cublas1.hpp>
#include <boost/numeric/bindings/cublas/cublas2.hpp>
#include <boost/numeric/bindings/cublas/cublas3.hpp>

#include <cublas.h>

namespace boost {
namespace numeric {
namespace bindings {
namespace cublas {

void shutdown() {
    cublasShutdown();
    check_status();
}

void init() {
    cublasInit();
    try {
	check_status();
    }
    catch (cublas_runtime_error e) {
	throw cublas_status_init_failed(e.what());
    }
}

}
}
}
}

#endif // boost_cublasbindings_cublasbindings_hpp
