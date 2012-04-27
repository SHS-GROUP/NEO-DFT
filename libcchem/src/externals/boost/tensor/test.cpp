#include "tensor/tensor.hpp"
#include "tensor/permutation.hpp"
#include "tensor/swap.hpp"
#include "tensor/permute.hpp"
#include "tensor/contract.hpp"
#include "tensor/operator.hpp"

// #include <boost/typeof/typeof.hpp>

// template<class A_, class B_, class C_>
// void function(tensor::expression<A_> &a_,
// 	      tensor::expression<B_> &b_,
// 	      const tensor::expression<C_> &c_) {

//     BOOST_AUTO(&A, a_());
//     BOOST_AUTO(&B, b_());
//     BOOST_AUTO(const &C, c_());

//     int i = 9;

//     tensor::index<0> a;
//     tensor::index<1> b;
//     tensor::index<2> c;
//     tensor::index<3> d;
//     tensor::index<4> e;

//     A(i,i,i) = C(b,c,a,i)*A(a,b,c);
//     A(i,i,i) = C(b,c,a,i)*(A(a,b,c) + B(a,b,i,c));
//     B(i,a,i,i) = A(a,b,i)*A(b,i,i);


// }

int main() {

    tensor::Tensor<3> A(10,20,10);
    size_t size[] = { 10, 10, 20, 10 };
    tensor::Tensor<4> B(size), C(size);

    //  for (int i = 0; i < A.size(0); ++i) {
    // 	A(i,i,i) = i;
    // }

    // 	     for (int i = 0; i < A.size(0); ++i) {
    // 		     std::cout << A(i,i,i)  << std::endl;
    // 		 }

    //int i = 0;

    tensor::index<'a'> a;
    tensor::index<'b'> b;
    tensor::index<'c'> c;
    tensor::index<'d'> d;
    tensor::index<'e'> e;
    tensor::index<'i'> i;
    tensor::index<'j'> j;

    C[0][0][0][0] = 6;

    BOOST_AUTO(I, tensor::Operator::I);
    BOOST_AUTO(P, tensor::Operator::P);

    (I - P(j,i))(I - P(i,j))(A(i,j,0));

    swap(A(b,0,0), A(0,b,0));
    tensor::permute<1,0>(A(b,0,c));
	 
    //B(a,0,0,0);
    C(a,c,0,d);

    // //     tensor::
    contract(1, A(b,0,a), B(a,c,d,0), 0, C(b,c,d,0));
    A(b,c,d) = contract(B(b,0,0,a), C(a,c,d,0));
    A(b,c,d) = 5.0*contract(B(b,a,0,0), C(a,c,d,0));

    //permutation(a,b,c)(A(a,b,c));

    // A(b,c,d) -= 5*contract(B(b,0,a,0), C(a,c,d,0));
    // A(b,c,d) += 4*contract(B(b,0,a,0), C(a,c,d,0));

}

//     int i = 9;


//     A(c,d,i) = C(b,c,a,i)*A(a,b,i);
//     A(i,i,i) = C(b,c,a,i)*(A(a,b,c) + B(a,b,i,c));
//     B(i,a,i,i) = A(a,b,i)*A(b,i,i);

//     function(A,B, A);

//     (A(a,b,c))=C(b,c,a,i);

//     (A(a,i,i) - A(a,i,i))(a);

//     (-A(a,b,d))(i,i,i);
//     //A(i,i,i);

//     //A(a,b,d) + B(a,b,i,d);
//     // T(b,c) = C(b,a)*D(a,c);


// }

