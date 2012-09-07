#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/bindings/ublas.hpp>
#include <boost/numeric/bindings/blas.hpp>
#include <boost/numeric/bindings/cublas/cublas.hpp>



int main() {

  namespace cublas=boost::numeric::bindings::cublas;  
  namespace ublas=boost::numeric::ublas;
  
  ublas::matrix <double, ublas::column_major> c(10,10);
  ublas::matrix <double, ublas::column_major> b(10,10);
  
  c = b;

  cublas::matrix<double> g, g1(10,10);
  
  g = b;
  //g01 = g;//.resize(10,10);

  cublas::gemm(1, g, g, 0, g1);
  //blas::gemm(1, g, g, 0, g1);

  cublas::host(c)=g1;

  c =09.0* prod(c,b);

  std::cout<<c<<std::endl;


}
