#ifndef _CXX_MATH_MATH_HPP_
#define _CXX_MATH_MATH_HPP_

#ifndef __host__
#define __host__
#endif

#ifndef __device__
#define __device__
#endif

namespace cxx {
    namespace math {

	template<typename T, size_t N>
	__host__ __device__
	T dot(const T (&v)[N]) {
	    T d = 0;
	    for (size_t i = 0; i < N; ++i) { d += v[i]*v[i]; }
	    return d;
	}

	template<typename T, typename U, size_t N>
	__host__ __device__
	T multiply(const U (&v)[N], const T&) {
	    T d = v[0];
	    for (size_t i = 1; i < N; ++i) { d *= v[i]; }
	    return d;
	}

	template<typename T, size_t N>
	__host__ __device__
	T multiply(const T (&v)[N]) {  return multiply(v, T()); }

	template<typename T, size_t N>
	__host__ __device__
	T multiply(const boost::array<T,N> &v) {  return multiply(v.elems, T()); }

	template<typename T, size_t N>
	__host__ __device__
	T add(const T (&v)[N]) {
	    T d = T(0);
	    for (size_t i = 0; i < N; ++i) { d += v[i]; }
	    return d;
	}


    }
}

#endif /* _CXX_MATH_MATH_HPP_ */
