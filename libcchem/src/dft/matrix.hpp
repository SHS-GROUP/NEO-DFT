#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <boost/typeof/typeof.hpp>
#include <boost/mpl/size_t.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/io.hpp>


#define BOOST_PROFILE_ENABLE
#include "boost/utility/profiler.hpp"


namespace block {
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
		typedef T type;
#else
		typedef T type __attribute__((aligned (8)));
#endif
	    };

	    template<class A, class B, class C>
	    static void
	    multiply(const A (&a)[L*L], const B (&b)[L*L], C (&c)[L*L]) {
		multiply_impl<A,B,C>(a, b, c);
	    }


	    template<class A, class B, class C>
	    static void 
	    multiply_impl(const typename aligned<A>::type* __restrict__ a,
			  const typename aligned<B>::type* __restrict__ b,
			  typename aligned<C>::type* __restrict__ c) {

		typedef typename aligned<A>::type A_;
		typedef typename aligned<B>::type B_;
		typedef typename aligned<C>::type C_;

		for (size_t j = 0; j < L; ++j) {
		    const B_* __restrict__ bj = b + j*L;
		    C_ c_[L] = { 0 };
		    for (size_t k = 0; k < L; ++k) {
			const A_* __restrict__ ak = a+k*L;
			B_ bjk = bj[k];
			for (size_t i = 0; i < L; ++i) {
			    c_[i] += bjk*ak[i];
			}
		    }
		    C_* __restrict__ cj = c + j*L;
		    for (size_t i = 0; i < L; ++i) {
			cj[i] += c_[i];
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
	    template<size_t N, size_t L, typename U, class A>
	    void operator()(const U &alpha, const A (&a)[L],
			    const ublas::range &r1,
			    const ublas::range &r2) const {
		block<N>::store(alpha, a, beta, ublas::project(c, r1, r2));
	    }
	};


    } // namespace detail
} // namespace block

namespace block {

    template<size_t L, typename T, class A, class B, class F>
    void prod(T alpha,
	      const boost::numeric::ublas::matrix_expression<A> &a,
	      const boost::numeric::ublas::matrix_expression<B> &b,
	      const F &f, double cutoff = 0) {

	BOOST_PROFILE_FUNCTION();

	namespace ublas = boost::numeric::ublas;
	typedef detail::block<L> block;
	typedef typename block::range range;

	int M = int(a().size1());
	int K = int(a().size2());
	int N = int(b().size2());

	typename block::max test(ublas::trans(a()), b(), cutoff/alpha);

#pragma omp parallel for schedule(dynamic,1)
	for (int j = 0; j < N; j += L) {
	    range rj(j, N);

	    for (int i = 0; i < M; i += L) {
		typename block::template
		    aligned<typename F::value_type>::type c_[L*L] = { 0 };
		range ri(i, M);
		for (int k = 0; k < K; k += L) {
		    range rk(k, K);

		    if (!test.ab(k,i,k,j)) continue;

		    typename block::template
			aligned<typename A::value_type>::type a_[L*L] = { 0 };
		    typename block::template
			aligned<typename B::value_type>::type b_[L*L] = { 0 };

		    block::load(a_, ublas::project(a(), ri, rk)); 
		    block::load(b_, ublas::project(b(), rk, rj));

		    block::multiply(a_, b_, c_);
		}
		f.template operator()<L,L*L>(alpha, c_, ri, rj);
	    }
	}

    }


    template<size_t L, typename T, class A, class B, typename U, class C>
    void prod(T alpha,
	      const boost::numeric::ublas::matrix_expression<A> &a,
	      const boost::numeric::ublas::matrix_expression<B> &b,
	      U beta,
	      boost::numeric::ublas::matrix_expression<C> &c,
	      double cutoff = 0) {
	prod<L>(alpha, a, b, detail::store_functor<U,C>(beta, c()), cutoff);
    }


    // C(i,j) = A(i,k)*B(j,k)' + B(i,k)*A(j,k)'
    template<size_t L, typename T, class A, class B, typename U, class C>
    void outer_prod(T alpha,
		    const boost::numeric::ublas::matrix_expression<A> &a,
		    const boost::numeric::ublas::matrix_expression<B> &b,
		    U beta,
		    boost::numeric::ublas::matrix_expression<C> &c,
		    double cutoff = 0) {

	BOOST_PROFILE_FUNCTION();

	namespace ublas = boost::numeric::ublas;
	typedef detail::block<L> block;
	typedef typename block::range range;

	int N = int(a().size1());
	int K = int(a().size2());

	typename block::max test(a(), b(), cutoff/alpha);

	assert((beta == 1.0) && "not implemented yet");

#pragma omp parallel for schedule(dynamic,1)
	for (int j_ = 0; j_ < N; j_ += L) {
	    int j = N - std::min<int>(j_+N%L,N); // long column first
	    for (int k = 0; k < K; k += L) {

		range rj(j, N);
		range rk(k, K);

		typename block::template
		    aligned<typename A::value_type>::type a_kj[L*L] = { 0 };
		block::load(a_kj, ublas::project(ublas::trans(a()), rk, rj));

		typename block::template
		    aligned<typename B::value_type>::type b_kj[L*L] = { 0 };
		block::load(b_kj, ublas::project(ublas::trans(b()), rk, rj));

		for (int i = 0; i <= j; i += L) {
		    range ri(i, j+rj.size());

		    typename block::template
			aligned<typename C::value_type>::type c_ij[L*L] = { 0 };

		    bool store = false;
		    if (test.ab(i,k,j,k) && test.a(i,k,b_kj)) {
			typename block::template
			    aligned<typename A::value_type>::type a_ik[L*L] = { 0 };
			block::load(a_ik, ublas::project(a(), ri, rk));
			block::multiply(a_ik, b_kj, c_ij);
			store = true;
		    }

		    if (test.ba(i,k,j,k) && test.b(i,k,a_kj)) {
			typename block::template
			    aligned<typename B::value_type>::type b_ik[L*L] = { 0 };
			block::load(b_ik, ublas::project(b(), ri, rk)); 
			block::multiply(b_ik, a_kj, c_ij);
			store = true;
		    }

		    if (!store) continue;
		    block::store(alpha, c_ij, 1.0, ublas::project(c(), ri, rj));

		}
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

}

#endif // MATRIX_HPP 
