#ifndef _CUDA_VECTOR_OPERATOR_HPP_
#define _CUDA_VECTOR_OPERATOR_HPP_

#include <iostream>
#include <string>

#include <cuda.h>
#include <cuda_runtime.h>

namespace cuda {

    template<class V>
    struct type_traits {};
	
    template<>
    struct type_traits< dim3> {
	typedef  unsigned int value_type;
	static const size_t size = 3;
    };

    template<>
    struct type_traits<int4> {
	typedef int value_type;
	static const size_t size = 4;
    };

    template<>
    struct type_traits<ushort3> {
	typedef ushort value_type;
	static const size_t size = 3;
    };

    template<class V>
    struct type_traits<const V> : type_traits<V> {};

    template<class Vector, size_t N>
    struct wrapper_base {
	typedef typename type_traits<Vector>::value_type value_type;
	Vector &v;
	__host__ __device__
	wrapper_base(Vector &vector) : v(vector) {}
	size_t size() const { return N ; }
	__host__ __device__
	value_type& operator[](size_t i) { return *(&v.x + i); }
	__host__ __device__
	const value_type& operator[](size_t i) const { return *(&v.x + i); }
    };

    template<class V, size_t N>
    struct wrapper;


    template<class Vector>
    struct wrapper<Vector, 2> :  wrapper_base<Vector, 2> {
	typedef wrapper_base<Vector, 2> base;
	typedef typename base::value_type value_type;
	__host__ __device__
	wrapper(Vector &vector) : base(vector) {}
	value_type multiply() const { return (base::v.x*base::v.y); }
    };

    template<class Vector>
    struct wrapper<Vector, 3> :  wrapper_base<Vector, 3> {
	typedef wrapper_base<Vector, 3> base;
	typedef typename base::value_type value_type;
	__host__ __device__
	wrapper(Vector &vector) : base(vector) {}
	value_type multiply() const { return (base::v.x*base::v.y*base::v.z); }
    };

    template<class Vector>
    struct wrapper<Vector, 4> :  wrapper_base<Vector, 4> {
	typedef wrapper_base<Vector, 4> base;
	typedef typename base::value_type value_type;
	__host__ __device__
	wrapper(Vector &vector) : base(vector) {}
	template<class A>
	__host__ __device__
	wrapper& operator=(const A &array) {
	    base::v.x = array[0];
	    base::v.y = array[1];
	    base::v.z = array[2];
	    base::v.w = array[3];
	    return *this;
	}
	template<typename T>
	__host__ __device__
	void __copy(T (&array)[4]) const {
	    array[0] = base::v.x;
	    array[1] = base::v.y;
	    array[2] = base::v.z;
	    array[3] = base::v.w;
	}
	value_type multiply() const { return (base::v.x*base::v.y*base::v.z*base::v.w); }
    };

    template<class V>
    __host__ __device__
    struct wrapper<V, type_traits<V>::size> wrap(V &v) {
	return wrapper<V, type_traits<V>::size>(v);
    }

    namespace detail {
	template<class V>
	// __host__ __device__
	std::ostream& operator<<(std::ostream &os, V vector) {
	    using std::operator<<;
	    cuda::wrapper<V, cuda::type_traits<V>::size> w(vector);
	    os << "{ ";
	    for (size_t i = 0; i < w.size()-1; ++i) os << w[i] << ", ";
	    return os << w[w.size()-1] << " }";
	}
    }

    template<class V>
    __host__ __device__
    void println(V &vector) {
	cuda::wrapper<V, cuda::type_traits<V>::size> w(vector);
#ifdef __DEVICE_EMULATION__
	printf("{ ");
	for (size_t i = 0; i < w.size()-1; ++i) printf("%i, ", w[i]);
        printf("%i }\n", w[w.size()-1]);
#endif
    }

    template<class V, typename T,size_t N>
    __host__ __device__
    void copy(const wrapper<V,N> &w, T (&array)[N]) {
	w.__copy(array);
    }

}


inline std::ostream& operator<<(std::ostream &os, const ushort3 &vector) {
    return cuda::detail::operator<<(os, vector);
}

inline std::ostream& operator<<(std::ostream &os, const int4 &vector) {
    return cuda::detail::operator<<(os, vector);
}

inline std::ostream& operator<<(std::ostream &os, const dim3 &vector) {
    return cuda::detail::operator<<(os, vector);
}

template<class V, typename T>
__host__ __device__
void copy(const V &vector, T (&array)[cuda::type_traits<V>::size]) {
    cuda::copy(cuda::wrap(vector), array);
}

template<class V>
typename cuda::type_traits<V>::value_type multiply(const V &vector) {
    return cuda::wrap(vector).multiply();
}

#endif /* _CUDA_VECTOR_OPERATOR_HPP_ */
