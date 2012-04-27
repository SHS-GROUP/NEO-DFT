#ifndef BOOST_CUBLAS_BINDINGS_CUBLAS_HPP
#define BOOST_CUBLAS_BINDINGS_CUBLAS_HPP

#include <boost/numeric/bindings/cublas/forward.hpp>
#include <boost/numeric/bindings/cublas/vector.hpp>
#include <boost/numeric/bindings/cublas/matrix.hpp>
#include <boost/numeric/bindings/cublas/host.hpp>
#include <boost/numeric/bindings/cublas/expression.hpp>
#include <boost/numeric/bindings/cublas/exception.hpp>

#include <boost/numeric/bindings/cublas/level1.hpp>
#include <boost/numeric/bindings/cublas/level2.hpp>
#include <boost/numeric/bindings/cublas/level3.hpp>

#include <cublas.h>

namespace boost {
namespace numeric {
namespace bindings {
namespace cublas {

inline void shutdown() {
    cublasShutdown();
    check_status();
}

inline void init() {
    cublasInit();
    try {
	check_status();
    }
    catch (cublas_runtime_error e) {
	throw cublas_status_init_failed(e.what());
    }
}

  inline void set_stream(cudaStream_t stream) {
#ifdef cublasSetStream
    cublasSetStream(stream);
#else
    cublasSetKernelStream(stream);
#endif
  }

}
}
}
}

#endif // BOOST_CUBLAS_BINDINGS_CUBLAS_HPP
