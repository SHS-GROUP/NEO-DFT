#ifndef RYSQ_KERNEL_VECTOR_HPP
#define RYSQ_KERNEL_VECTOR_HPP

#include <boost/array.hpp>
#include <boost/config.hpp>

namespace rysq {
namespace kernel {

    template<size_t N>
    struct vector_operator {
	template<class V, class A>
	BOOST_GPU_ENABLED
	static void assign(V &v, const A &a) {
	    for (size_t i = 0; i < N; ++i) {
		v[i] = a[i];
	    }
	}
	template<class V, class A>
	BOOST_GPU_ENABLED
	static void multiply_assign(V &v, const A &a) {
	    for (size_t i = 0; i < N; ++i) {
		v[i] *= a;
	    }
	}
	template<class V, class A>
	BOOST_GPU_ENABLED
	static void divide_assign(V &v, const A &a) {
	    for (size_t i = 0; i < N; ++i) {
		v[i] /= a;
	    }
	}
	template<class V, class A>
	BOOST_GPU_ENABLED
	static void minus_assign(V &v, const A &a) {
	    for (size_t i = 0; i < N; ++i) {
		v[i] -= a[i];
	    }
	}
    };

    template<size_t N, typename T, class A, class V>
    struct vector_base {
	BOOST_GPU_ENABLED
	T* data() { return elems; }
	BOOST_GPU_ENABLED
	const T* data() const { return elems; }
	BOOST_GPU_ENABLED
	T& operator[](size_t i)  { return elems[i]; }
	BOOST_GPU_ENABLED
	const T& operator[](size_t i) const { return elems[i]; }
	BOOST_GPU_ENABLED
	V operator-(const V &rhs) const {
	    V v;
	    for (size_t i = 0; i < N; ++i) {
		v[i] = elems[i] - rhs[i];
	    }
	    return v;
	}
	BOOST_GPU_ENABLED
	operator A&() { return elems; }
	BOOST_GPU_ENABLED
	operator const A&() const { return elems; }

	BOOST_GPU_ENABLED
	T dot() const {
	    T q = 0.0;
	    for (size_t i = 0; i < N; ++i) {
		q += elems[i]*elems[i];
	    }
	    return q;
	}

	BOOST_GPU_ENABLED
	static V center(T w1, const V &v1, T w2, const V &v2) {
	    V v;
	    T w = T(1)/(w1 + w2);
	    for (size_t i = 0; i < N; ++i) {
		v[i] = (w1*v1[i] + w2*v2[i])*w;
	    }
	    return v;
	}
	A elems __attribute__ ((aligned(16)));

    protected:
	template<class C>
	BOOST_GPU_ENABLED
	void assign(const C &v) {
	    for (size_t i = 0; i < N; ++i) {
		elems[i] = v[i];
	    }
	}
    };

    template<size_t N, typename T = double>
    struct vector : vector_base<N, T, T[N], vector<N,T> > {
    private:
	typedef vector_base<N, T, T[N], vector<N,T> > base;
    public:

	struct adapter : vector_base<N, T, T*, vector> {
	    typedef vector_base<N, T, T*, vector> base;
	    BOOST_GPU_ENABLED
	    explicit adapter(T *data) { base::elems = data; }
	    using base::operator-;
	    BOOST_GPU_ENABLED
	    vector operator-(const adapter &rhs) const {
		vector v(this->elems);
		vector_operator<N>::minus_assign(v, rhs);
		return v;
	    }
	};

	BOOST_GPU_ENABLED
	vector() {} 

	BOOST_GPU_ENABLED
	void clear() {
	    for (int i = 0; i < N; ++i) {
		(*this)[i] = 0;
	    }
	}

	BOOST_GPU_ENABLED
	vector(const adapter &array) {
	    base::assign(array);
	}

	BOOST_GPU_ENABLED
	explicit vector(const boost::array<T,N> &array) {
	    base::assign(array.elems);
	}

	BOOST_GPU_ENABLED
	explicit vector(const T *array) {
	    base::assign(array);
	}

	BOOST_GPU_ENABLED
	vector& operator*=(const T &v) {
	    vector_operator<N>::multiply_assign(*this, v);
	    return *this;
	}

	BOOST_GPU_ENABLED
	vector& operator/=(const T &v) {
	    vector_operator<N>::divide_assign(*this, v);
	    return *this;
	}

    };

}
}

#endif /* RYSQ_KERNEL_VECTOR_HPP */
