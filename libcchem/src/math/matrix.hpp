#ifndef MATH_MATRIX_HPP
#define MATH_MATRIX_HPP

#include <boost/typeof/typeof.hpp>
#include <boost/mpl/size_t.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "boost/utility/profiler.hpp"

#include "Eigen/Dense"

namespace math {
namespace matrix {
namespace detail {

    namespace ublas = boost::numeric::ublas;


    template<class E>
    struct iterator {
	typedef typename E::iterator2 column;
	typedef typename E::iterator1 row;
    };

    template<class E>
    struct iterator<const E> {
	typedef typename E::const_iterator2 column;
	typedef typename E::const_iterator1 row;
    };


    template<size_t L>
    struct block {
    private:

	struct load_impl {
	    template<typename T, typename U>
	    void operator()(T &t, const U &u) const  { t = u; }
	};

	struct store_impl1 {
	    template<typename T, typename U>
	    void operator()(const T &t, U &u) const  { u = tx*t; }
	    explicit store_impl1(double tx) : tx(tx) {}
	    double tx;
	};
 
	struct store_impl2 {
	    template<typename T, typename U>
	    void operator()(const T &t, U &u) const  {
		u = ux*u + tx*t;
	    }
	    store_impl2(double tx, double ux) : tx(tx), ux(ux) {}
	    double tx, ux;
	};

	template<class F, class A, class E>
	static void apply(const F &f, A &a, E &e) {
	    typedef typename iterator<E>::column column;
	    typedef typename iterator<E>::row row;
	    BOOST_AUTO(aj, &a[0]);
	    for (column c = e.begin2(); c != e.end2(); ++c) {
		int i = 0;
		for (row r = c.begin(); r != c.end(); ++r) {
		    f(aj[i++], *r);
		}
		aj += L;
	    }
	}

    public:

	template<class A, class B>
	static void load(A (&a)[L][L], const B b) {
	    A *a_ = &a[0][0];
	    load(a_, b);
	}

	template<class A, class B>
	static void load(A &a, const B b) {
	    apply(load_impl(), a, b);
	}

	template<class A, class B>
	static void store(double ax, const A &a, double bx, B b) {
	    if (typename B::value_type(1)*bx == typename B::value_type(0)) {
		apply(store_impl1(ax), a, b);
	    }
	    else {
		apply(store_impl2(ax, bx), a, b);
	    }
	}

	template<class A, class B>
	static void store(double ax, const A (&a)[L][L], double bx, B b) {
	    store(ax, &a[0][0], bx, b);
	}

	struct range : ublas::range {
	    range(size_t i, size_t n)
		: ublas::range(i, i + std::min(L, n-i)) {}
	};

	struct max {
	    typedef ublas::matrix<double, ublas::column_major> matrix;
	    template<class A, class B>
	    max(const A &a, const B &b, double value) {
		initialize(a, a_);
		initialize(b, b_);
		value_ = value/L;
	    }
	    bool ab(size_t i, size_t j, size_t k) const {
		return (a(i,k)*b(k,j) > value_);
	    }
	    double ba(size_t i, size_t j, size_t k) const {
		return (b(i,k)*a(k,j) > value_);
	    }
	    bool ab(size_t i, size_t j, size_t k, size_t l) const {
		return (a(i,j)*b(k,l) > value_);
	    }
	    bool ba(size_t i, size_t j, size_t k, size_t l) const {
		return (a(k,l)*b(i,j) > value_);
	    }
	    double a(size_t i, size_t j) const { return a_(i/L, j/L); }
	    double b(size_t i, size_t j) const { return b_(i/L, j/L); }

	    template<class A>
	    bool a(size_t i, size_t j, const A (&m)[L*L]) const {
		return test(a(i,j), m, value_);
	    }
	    template<class A>
	    bool b(size_t i, size_t j, const A (&m)[L*L]) const {
		return test(b(i,j), m, value_);
	    }

	private:
	    matrix a_, b_;
	    double value_;

	    template<class A>
	    static void initialize(const A &a, matrix &norm) {
		size_t M = a.size1();
		size_t N = a.size2();
		norm.resize(M/L + 1, N/L + 1, false);
		for (size_t j = 0; j < N; j += L) {
		    range rj(j, N);
		    for (size_t i = 0; i < M; i += L) {
			range ri(i, M);
			norm(i/L, j/L) =
			    ublas::norm_inf(ublas::project(a, ri, rj));
		    }
		}
	    }

	    template<typename T, class A, typename U>
	    static bool test(T v, const A (&m)[L*L], U cutoff) {
		A w[L] = { 0 };
		for (size_t i = 0; i < L; ++i) {
		    for (size_t j = 0; j < L; ++j) {
			w[j] += v*m[j+i*L];
		    }
		}
		int n = 0;
		for (size_t j = 0; j < L; ++j) {
		    n += (std::abs(w[j]) > cutoff);
		}
		return (n > 0);
	    }

	};

	template<class A>
	static void transpose(A &a) {
	    for (size_t j = 0; j < L; ++j) {
		for (size_t i = 0; i < j; ++i) {
		    a[j][i] = a[i][j];
		}
	    }
	}

	template<class T>
	struct aligned {
#if (defined(__GNUC__) && ((__GNUC__ < 4) || (__GNUG__ == 4 && __GNUC_MINOR__ < 3)))
#warning "alignment not implemented for GNUC < 4.3"
	    typedef T type[L][L];
#else
	    typedef T type[L][L] __attribute__((aligned (16)));
#endif
	    typedef T value_type;
	};

	template<class A, class B, class C>
	static void
	multiply(const A &a, const B &b, C &c) {
		 // typename aligned<A>::type &a,
		 // const typename aligned<B>::type &b,
		 // typename aligned<C>::type &c) {
	    //const A (&a)[L][L], const B (&b)[L][L], C (&c)[L][L]) {
	    multiply_impl<double,double,double>(a, b, c);
	}


	template<class A, class B, class C>
	static void 
	multiply_impl(const typename aligned<A>::type &a,
		      const typename aligned<B>::type &b,
		      typename aligned<C>::type &c) {

	    typedef typename aligned<A>::value_type A_;
	    typedef typename aligned<B>::value_type B_;
	    typedef typename aligned<C>::value_type C_;

	    // for (size_t i = 0; i < L; ++i) {
	    // 	for (size_t j = 0; j < L; ++j) {
	    // 	    C_ c_ = 0;
	    // 	    const A_* __restrict__ a_ = a+i*L;
	    // 	    const B_* __restrict__ b_ = b+j*L;
	    // 	    for (size_t k = 0; k < L; ++k) {
	    // 		c_ += b_[k]*a_[k];
	    // 	    }
	    // 	    c[i+j*L] += c_;
	    // 	}
	    // }
	    // return;

	    for (size_t j = 0; j < L; ++j) {
		C_ c_[L] __attribute__((aligned (16))) = { 0 } ;
		for (size_t k = 0; k < L; k += 2) {
		    B_ bjk[2] __attribute__((aligned (16))) =
			{ b[j][k+0], b[j][k+1] };
		    if (fabs(bjk[0])+fabs(bjk[1]) < 1e-10) continue;
		    for (size_t i = 0; i < L; i += 1) {
			c_[i+0] += bjk[0]*a[k+0][i];
			c_[i+0] += bjk[1]*a[k+1][i];
		    }
		}
		for (size_t i = 0; i < L; ++i) {
		    c[j][i] += c_[i];
		}
	    }

	}

    };

    template<typename T, class C>
    struct store_functor {
	typedef typename C::value_type value_type;
	T beta; C &c;
	store_functor(T beta, C &c) : beta(beta), c(c) {}  
	template<size_t N>
	struct squared { static const size_t value = N*N; };
	template<size_t N, typename U, class A>
	void operator()(const U &alpha, const A (&a)[N][N],
			const ublas::range &r1,
			const ublas::range &r2) const {
	    block<N>::store(alpha, a, beta, ublas::project(c, r1, r2));
	}
    };


    template<size_t L>
    struct block_max {
	typedef ublas::matrix<float, ublas::column_major> matrix;
	template<class A>
	explicit block_max(const boost::numeric::ublas::matrix_expression<A> &a) {
	    size_t m = a().size1();
	    size_t n = a().size2();
	    a_.resize((m+L-1)/L, (n+L-1)/L, false);
	    for (size_t j = 0; j < n; j += L) {
		ublas::range rj(j, std::min(j+L,n));
		for (size_t i = 0; i < m; i += L) {
		    ublas::range ri(i, std::min(i+L,m));
		    a_(i/L, j/L) = ublas::norm_inf(ublas::project(a(), ri, rj));
		}
	    }
	}
	float operator()(int i, int j) const {
	    return a_(i/L,j/L);
	}
    private:
	matrix a_;
    };



    struct load {
	template<class A, class B>
	void operator()(A &a, const B &b) const { a = b; }
    };

    template<typename T = void, typename U = void>
    struct store {
	T alpha;
	U beta;
	explicit store(const T &alpha, const U &beta)
	    : alpha(alpha), beta(beta) {}
	template<class A, class B>
	void operator()(const A &a, B &b) const {
	    b = beta*b + alpha*a;
	    //b = a;
	}
    };

    template<>
    struct store<> {
	template<class A, class B>
	void operator()(const A &a, B &b) const {
	    b = a;
	}
    };

    template<typename T>
    struct scalar {
	T value;
	explicit scalar(const T &value) : value(value) {}
	const T& operator()(int i = 0, int j = 0) const { return value; }
    };

    template<typename T, size_t L>
    struct vector : Eigen::Matrix<T,L,1> {
	vector() { clear(); }
	void clear() {
	    for (int i = 0; i < L; ++i) {
		(*this)[i] = 0;
	    }
	}
    };

    template<typename T, size_t L>
    struct matrix : Eigen::Matrix<T,L,L> {
	typedef typename Eigen::MatrixBase<
	    Eigen::Matrix<T,L,L> >::ColXpr column;
	typedef typename Eigen::MatrixBase<
	    const Eigen::Matrix<T,L,L> >::ColXpr const_column;
	matrix() { clear(); }
	void clear() {
	    for (int j = 0; j < L; ++j) {
		for (int i = 0; i < L; ++i) {
		    (*this)(i,j) = 0;
		}
	    }
	}
	column operator[](int i) { return this->col(i); }
	const_column operator[](int i) const { return this->col(i); }
    };

    // template<typename T, size_t L>
    // struct matrix {
    // 	T data_[L][L];
    // 	matrix() { clear(); }
    // 	T& operator()(int i, int j) { return data_[j][i]; }
    // 	const T& operator()(int i, int j) const { return data_[j][i]; }
    // 	void clear() {
    // 	    for (int j = 0; j < L; ++j) {
    // 		for (int i = 0; i < L; ++i) {
    // 		    data_[j][i] = 0;
    // 		}
    // 	    }
    // 	}
    // };

    template<class F, class A, class B>
    inline void apply(const F &f, A &a, B &b) {
	typedef BOOST_TYPEOF(b.begin2()) column;
	int j = 0;
	for (column c = b.begin2(); c < b.end2(); ++c, ++j) {
	    typedef BOOST_TYPEOF(c.begin()) row; 
	    int i = 0;
	    for (row r = c.begin(); r < c.end(); ++r, ++i) {
		f(a(i,j), *r);
	    }
	}
    }

    template<size_t L, class A, class B, class C>
    inline void multiply(const matrix<A,L> &a,
			 const matrix<B,L> &b,
			 matrix<C,L> &c,
			 B cutoff) {
	for (size_t j = 0; j < L; ++j) {
	    // vector<C,L> c_;
	    for (size_t k = 0; k < L; k += 1) {
		// B bjk[2] = { b(k+0,j), b(k+1,j) };
		// if (fabs(bjk[0])+fabs(bjk[1]) < cutoff) continue;
		// {
		//     c_ += bjk[0]*a[k+0];
		//     c_ += bjk[1]*a[k+1];
		// }
		if (fabs(b(k+0,j)) < cutoff) continue;
		c[j] += b(k+0,j)*a[k+0];
		// for (int i = 0; i < L; ++i) {
		//     c_[i] += bjk[0]*a[k+0][i];
		//     c_[i] += bjk[1]*a[k+1][i];
		// }
	    }
	    // for (size_t i = 0; i < L; ++i) {
	    // 	c[j][i] += c_[i];
	    // }
	    // c[j] += c_;
	}
    }


} // namespace detail
} // namespace matrix
} // namespace math


namespace math {
namespace matrix {

    template<int L, typename T, class A, class B, typename U, class C>
    void prod(T alpha,
	      const boost::numeric::ublas::matrix_expression<A> &a,
	      const boost::numeric::ublas::matrix_expression<B> &b,
	      U beta,
	      boost::numeric::ublas::matrix_expression<C> &c,
	      double cutoff = 0) {

	BOOST_PROFILE_LINE;

	namespace ublas = boost::numeric::ublas;
	typedef detail::block<L> block;
	typedef typename block::range range;

	using std::min;

	int M = int(a().size1());
	int K = int(a().size2());
	int N = int(b().size2());

	// std::cout << M << ":" << ":" << N << ":" << K
	// 	  << " cutoff = " << cutoff << std::endl;

	int tasks = 0;
#pragma omp parallel
	for (int j = 0, ij = 0, task = -1; j < N; j += L) {
	    for (int i = 0; i < M; i += L, ++ij) {

#pragma omp critical
		if (task < 0) task = tasks++; // initialize

		if (ij != task) continue;
#pragma omp critical
		task = tasks++;

		detail::matrix<double, L> c_;
		range rj(j, min(j+L,N));
		range ri(i, min(i+L,M));

		for (int k = 0; k < K; k += L) {
		    range rk(k, min(k+L,K));

		    detail::matrix<double, L> a_, b_;
		    detail::load load;
		    detail::apply(load, a_, ublas::project(a(), ri, rk));
		    detail::apply(load, b_, ublas::project(b(), rk, rj));
		    detail::multiply(a_, b_, c_, cutoff);
		}

		BOOST_AUTO(cij, ublas::project(c(), ri, rj));

		if (beta == U(0)) {
		    detail::store<> store;
		    detail::scalar<U> scalar(0);
		    detail::apply(store, scalar, cij);
		}

		detail::store<T,U> store(alpha, beta);
		detail::apply(store, c_, cij);

	    }
	}
    }


    template<typename T, class A, class B, typename U, class C>
    void prod(T alpha,
	      const boost::numeric::ublas::matrix_expression<A> &a,
	      const boost::numeric::ublas::matrix_expression<B> &b,
	      U beta,
	      boost::numeric::ublas::matrix_expression<C> &c,
	      double cutoff = 0) {
	prod<2*16>(alpha, a, b, beta, c, cutoff);
    }



    // C(i,j) = A(i,k)*B(j,k)' + B(i,k)*A(j,k)'
    template<int L, typename T, class A, class B, typename U, class C>
    void outer_prod(T alpha,
		    const boost::numeric::ublas::matrix_expression<A> &a,
		    const boost::numeric::ublas::matrix_expression<B> &b,
		    U beta,
		    boost::numeric::ublas::matrix_expression<C> &c,
		    double cutoff = 0) {

	BOOST_PROFILE_LINE;

	using std::min;
	namespace ublas = boost::numeric::ublas;
	typedef detail::block<L> block;
	typedef typename block::range range;

	int N = int(a().size1());
	int K = int(a().size2());

	assert((beta == 1.0) && "not implemented yet");

	int tasks = 0;
#pragma omp parallel
	for (int j = 0, ij = 0, task = -1; j < N; j += L) {
	    for (int i = 0; i <= j; i += L, ++ij) {

#pragma omp critical
		if (task < 0) task = tasks++; // initialize

		if (ij != task) continue;
#pragma omp critical
		task = tasks++;

		range rj(j, std::min(j+L,N));
		range ri(i, std::min<int>(i+L, j+rj.size())); 
		
		detail::matrix<double, L> c_ij;

		for (int k = 0; k < K; k += L) {
 
		    range rk(k, std::min(k+L,K));
 
		    detail::matrix<double, L> a_kj, b_kj, a_ik, b_ik;

		    detail::load load;
		    detail::apply(load, a_ik, ublas::project(a(), ri, rk));
		    detail::apply(load, b_ik, ublas::project(b(), ri, rk));
		    detail::apply(load, a_kj, ublas::project(ublas::trans(a()), rk, rj));
		    detail::apply(load, b_kj, ublas::project(ublas::trans(b()), rk, rj));

		    detail::multiply(b_ik, a_kj, c_ij, cutoff);
		    detail::multiply(a_ik, b_kj, c_ij, cutoff);

		}

		BOOST_AUTO(c_, ublas::project(c(), ri, rj));
		detail::store<T,U> store(alpha, beta);
		detail::apply(store, c_ij, c_);

	    }
	}
    
    }


    // C(i,j) = A(i,k)*B(j,k)' + B(i,k)*A(j,k)'
    template<typename T, class A, class B, typename U, class C>
    void outer_prod(T alpha,
		    const boost::numeric::ublas::matrix_expression<A> &a,
		    const boost::numeric::ublas::matrix_expression<B> &b,
		    U beta,
		    boost::numeric::ublas::matrix_expression<C> &c,
		    double cutoff = 0) {


	outer_prod<16*2>(alpha, a, b, beta, c, cutoff);
	return;


	BOOST_PROFILE_LINE;

    	namespace ublas = boost::numeric::ublas;

    	int N = a().size1();
    	int K = a().size2();

    	// C(i,j) = A(i,k)*B(j,k)' + B(i,k)*A(j,k)'

    	struct {
    	    double a, b;
    	} max = { 0, 0 };

    	for (int j = 0; j < N; j += 1) {
    	    for (int k = 0; k < K; k += 1) {
    		ublas::range ri(0,j+1);

    		double a_kj = alpha*a()(j,k);
    		double b_kj = alpha*b()(j,k);

    		max.a = std::max(max.a, std::abs(a()(j,k)));
    		max.b = std::max(max.b, std::abs(b()(j,k)));
		
    		bool ab = std::abs(alpha*max.a*b_kj);
    		bool ba = std::abs(alpha*max.b*a_kj);

    		BOOST_AUTO(cj_, ublas::column(c(),j));
    		BOOST_AUTO(cj, ublas::noalias(ublas::project(cj_, ri)));
    		BOOST_AUTO(ak, ublas::project(ublas::column(a(),k), ri));
    		BOOST_AUTO(bk, ublas::project(ublas::column(b(),k), ri));
		
    		if ((ab + ba) > cutoff) cj += (b_kj*ak + a_kj*bk);
    		else if (ab > cutoff) cj += b_kj*ak;
    		else if (ba > cutoff) cj += a_kj*bk;

    	    }
    	}
    }

} // namespace matrix
} // namespace math

#endif // MATH_MATRIX_HPP 
