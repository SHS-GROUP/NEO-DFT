#ifndef BINDINGS_BINDINGS_HPP
#define BINDINGS_BINDINGS_HPP

#include "utility/any.hpp"
// #include "runtime.hpp"
// #include "array/hdf5.hpp"

namespace bindings {

    utility::any<int>& any();

    // template<typename T>
    // void new_array(size_t N, const T* dims, const char *name) {
    // 	runtime::rt().put(name, new Array<double>(N, dims, name, Array<>::HDF5()));
    // }

    // template<typename T>
    // void array_put(const double *buffer, const T *start, const T *stop,
    // 		   const char *name) {
    // 	runtime::rt().get<Array<double>*>(name)->put(buffer, start, stop);
    // }

    // template<typename T>
    // void array_get(double *buffer, const T *start, const T *stop,
    // 		   const char *name) {
    // 	runtime::rt().get<Array<double>*>(name)->get(buffer, start, stop);
    // }

    // template<typename T>
    // double array_at(const T *index, const char *name) {
    // 	return runtime::rt().get<Array<double>*>(name)->at(index);
    // }

}


#endif // BINDINGS_BINDINGS_HPP
