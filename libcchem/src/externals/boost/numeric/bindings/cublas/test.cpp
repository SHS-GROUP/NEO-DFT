#include <boost/numeric/bindings/cublas/matrix.hpp>
#include <boost/numeric/bindings/cublas/cublas.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

int main() {

    namespace cublas = boost::numeric::bindings::cublas;
    namespace ublas = boost::numeric::ublas;

    cublas::init();

    cublas::matrix<double> A;
    cublas::matrix<double> B;
    B = cublas::matrix<double>(5,5);
    A = B;

    ublas::matrix<double> h(cublas::host(B));
    ublas::matrix<double, ublas::column_major> f(5,5);
    cublas::host(f) = A;

    cublas::matrix<double> C = cublas::gemm(A, cublas::trans(B));
    B += cublas::gemm(A, C);
    cublas::gemm(3.0, B, C, 6.0, A);

    cublas::row(A,0) = cublas::row(B,1);
    cublas::row(A,0) = cublas::column(B,0);
    
    using ublas::range;
    cublas::submatrix(A, range(0,2), range(0,2)) = 
    	cublas::submatrix(B, range(0,2), range(0,2));

    range r1(0,5), r2(0,5);
    cublas::submatrix(A, r1, r2) = cublas::gemm(A, C);
    cublas::submatrix(A, r1, r2) += cublas::gemm(A, C);

    cublas::shutdown();

}
