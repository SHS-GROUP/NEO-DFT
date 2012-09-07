#ifndef EIGEN_JACOBY_HPP
#define EIGEN_JACOBY_HPP

#include <cmath>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>


#include <boost/typeof/typeof.hpp>

namespace eigen {

    struct jacoby {

	struct rotation {
	    rotation(double t, double c) : s(t*c), tau(s/(1.0+c)) {}
	    template<typename T>
	    void operator()(T &p, T &q) const {
		T g = p, h = q;
		p = p - s*(h + p*tau);
		q = q + s*(g - q*tau);
	    }
	    double s, tau;
	};

	template<class M0, class M1>
	static void rotate(M0 &A,  M1 &V, double t, int p, int q) {
	    rotation rotate(t, 1.0/sqrt(1+t*t));
	    int N = A.size1();
	    for (int j = 0; j < p; j++) rotate(A(j,p), A(j,q));
	    for (int j = p+1; j < q; j++) rotate(A(p,j), A(j,q));
	    for (int j = q+1; j < N; j++) rotate(A(p,j), A(q,j));
	    for (int j = 0; j < N; j++) rotate(V(j,p), V(j,q));
	    A(p,q) = 0;
	}


	static double t_value(double pq, double g, double h) {
	    using std::abs;
	    double t;
	    if ((g + abs(h)) == abs(h)) {
		t = pq/h;
	    }
	    else {
		double theta = 0.5*h/pq;
		t = 1.0/(abs(theta) + sqrt(1.0 + theta*theta));
		if (theta < 0.0) t = -t;
	    }
	    return t;
	}


	template<class M0, class V0, class M1>
	struct sweep {
	    sweep(M0 &A, V0 &w, M1 &V) : A(A), w(w), V(V)
	    {
		initialize();
	    }
	    void initialize() {
		size_t N = A.size1();
		z.resize(N, false);
		b.resize(N, false);
		V.clear();
		z.clear();
		for (size_t i = 0; i < N; ++i) {
		    b(i) = w(i) = A(i,i);
		    V(i,i) = 1;
		}		
	    }

	    template<class R>
	    double operator()(const R &range);
	private:
	    M0 &A; V0 &w; M1 &V;
	    boost::numeric::ublas::vector<double> z, b;
	    size_t iter;
	};

	template<bool Tri>
	struct range_base {
	    range_base(size_t start1, size_t stop1,
		       size_t start2, size_t stop2)
		: start1_(start1), stop1_(stop1), 
		  start2_(start2), stop2_(stop2) {}
	    size_t start() const { return start2_; }
	    size_t stop() const { return stop2_; }
	    size_t start(size_t j) const { return (Tri) ? (j+1): start1_; }
	    size_t stop(size_t j) const { return stop1_; }
	private:
	    size_t start1_, stop1_;
	    size_t start2_, stop2_;
	};

	typedef range_base<false> block_range;

	struct triangular_range : range_base<true> {
	    triangular_range(size_t N)
		: range_base<true>(0, N, 0, N) {}
	};

	template<class M0, class V0, class M1>
	static void apply(M0 &A, V0 &w, M1 &V) {
	    size_t N = A.size1();
	    sweep<M0,V0,M1> sweep(A, w, V);
	    for (int i = 0; i < 50; ++i) {
		double sum = sweep(triangular_range(N));
		if (sum == 0) return;
	    }
	    throw std::runtime_error("too many iterations");
	}

    };

    // template<class M0>
    // void block_diagonalize(M0 &A, size_t N) {
    // 	namespace ublas = boost::numeric::ublas;
    // 	size_t M = A.size1();
    // 	ublas::vector<double> w(M);
    // 	ublas::matrix<double> V(M,M);
    // 	jacoby::sweep<M0,ublas::vector<double>,ublas::matrix<double> > sweep(A, w, V);
	
    // 	// for (int i = 0; i < 2; ++i) {
    // 	    sweep(jacoby::triangular_range(M));
    // 	// }

    // 	double off = 0;
    // 	for (int i = 0; i < 14; ++i) {
    // 	    //sweep(jacoby::triangular_range(N));
    // 	    double sum = sweep(jacoby::block_range(N, M, 0, N));
    // 	    if (sum == 0) return;
    // 	    double r = 1.6;
    // 	    double test = (sum > 1) ? (pow(sum, r) < off) : (pow(off, r) > sum);
    // 	    if ((off > 0 && !test)) {
    // 		std::cout << "triangular: " << sum << " " << off << std::endl;
    // 		sweep(jacoby::triangular_range(M));
    // 		off = 0;
    // 	    }
    // 	else off = sum;
    // 	}
    // 	throw std::runtime_error("too many iterations");
	
    // 	// for (size_t j = 0; j < N; ++j) {
    // 	//     for (size_t i = 0; i < N; ++i) {
    // 	// 	jacoby::rotate(
    // 	// }
    // }

    // struct block_jacoby {


    // 	template<class M0, class V0, class M1>
    // 	static void apply(M0 &A, V0 &w, M1 &V) {
    // 	    boost::numeric::ublas::symmetric_matrix<double> A_ = A;
    // 	    boost::numeric::ublas::vector<double> w_(w);
    // 	    apply_(A_, w_, V);
    // 	    std::cout << w_ << std::endl;
    // 	    jacoby::apply(A, w, V);
    // 	    std::cout << w << std::endl;
    // 	    if ( A.size1() > 22) throw;
    // 	}

    // 	template<class M0, class V0, class M1>
    // 	static void apply_(M0 &A, V0 &w, M1 &V) {
    // 	    namespace ublas = boost::numeric::ublas;
    // 	    size_t N = A.size1();
    // 	    if (N > 22) {
    // 		std::cout << "divide" << std::endl;
    // 		size_t N1 = N/2, N2 = N - N/2;
    // 		ublas::symmetric_matrix<double> A1(N1), A2(N2);
    // 		ublas::vector<double> w1(N1), w2(N2);
    // 		ublas::matrix<double> V1(N1,N1), V2(N2,N2);
    // 		block_diagonalize(A, N/2);
    // 		A1 = ublas::subrange(A, 0, N1, 0, N1);
    // 		A2 = ublas::subrange(A, N1, N, N1, N);
    // 		apply_(A1, w1, V1);
    // 		apply_(A2, w2, V2);
    // 		ublas::subrange(w, 0, N1).assign(w1);
    // 		ublas::subrange(w, N1, N).assign(w2);
    // 	    }
    // 	    else {
    // 		jacoby::apply(A, w, V);
    // 	    }
    // 	}

    // }; 


    template<class M0, class V0, class M1> template<class R>
    double jacoby::sweep<M0,V0,M1>::operator()(const R &range) {

	using std::abs;
	const int N = A.size1();

	double sum = 0;		
	for (size_t j = range.start(); j < range.stop(); ++j) {
	    for (size_t i = range.start(j); i < range.stop(j); ++i) {
		sum += abs(A(j,i));
	    }
	}
	// std::cout << sum << std::endl;
	if (sum == 0.0) return 0.0;

	double threshold = (iter < 4) ? (0.2*sum/(N*N)) : 0;

	for (size_t p = range.start(); p < range.stop(); ++p) {
	    for (size_t q = range.start(p); q < range.stop(p); ++q) {
		// for (int p = 0; p < N-1; ++p) {
		//     for (int q = p+1; q < N; ++q) {
		//std::cout << p << " "  << q  << " "<< A(p,q) << std::endl;
		double g = 100.0*abs(A(p,q));
		if (iter > 4 && (abs(w(p))+g) == abs(w(p)) && 
		    (abs(w(q))+g) == abs(w(q))){
		    A(p,q) = 0;
		    continue;
		}
		if (abs(A(p,q)) <= threshold) continue;

		double t = t_value(A(p,q), g, w(q) - w(p));
		double h = t*A(p,q);
		rotate(A, V, t, p, q);
		z(p) -= h;
		z(q) += h;
		w(p) -= h;
		w(q) += h;
	    }
	}

	for (int p = 0; p < N; ++p) {
	    b(p) += z(p);
	    w(p) = b(p);
	    z(p) = 0.0;
	}

	++iter;
	return sum;
    }

}

#endif // EIGEN_JACOBY_HPP
