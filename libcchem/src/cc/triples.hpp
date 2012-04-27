// #ifndef CC_TRIPLES_HPP
// #define CC_TRIPLES_HPP

// #include "tensor/order.hpp"
// #include "tensor/expression.hpp"
// #include "tensor/index.hpp"
// #include "tensor/tensor.hpp"
// #include "tensor/product.hpp"

// #include <boost/array.hpp>
// #include <boost/typeof/typeof.hpp>

// namespace cc {


//     template<class E0, class E1, class E2>
//     void triples(int i, int j, int k,
// 		 const tensor::expression<E0> &t2_, 
// 		 const tensor::expression<E1> &VVVO_,
// 		 const tensor::expression<E2> &OVOO_,
// 		 tensor::tensor<3> &X) {

// 	//int N = 0;

// 	using tensor::index_names::a;
// 	using tensor::index_names::b;
// 	using tensor::index_names::c;
// 	using tensor::index_names::e;
// 	using tensor::index_names::c;
// 	using tensor::index_names::c;



// 	const E0 &t2 = t2_();
// 	const E1 &VVVO = VVVO_();
// 	//const E2 &OVOO = OVOO_();

// 	// boost::mpl::print<BOOST_TYPEOF(t2(a,e,i,j))>::type;
// 	// boost::mpl::print<BOOST_TYPEOF(VVVO(b,c,e,k))>::type;
// 	// X(a,b,c) += t2(a,c,b,j);
// 	// (t2(a,e,i,j)*VVVO(b,c,e,k));
// 	1*t2(a,e,i,j);
// 	X(a,c,b) += 2.0*(t2(a,e,i,j)*VVVO(b,c,e,k));
// 	//	boost::mpl::print<BOOST_TYPEOF(X(a,b,c))>::type;

// 	// X(a,b,c) += (t2(a,e,i,j)*VVVO(b,c,e,k));

// 	// for (int c = 0; c < N; ++c) {
// 	//     X(a,b,c) += t2(a,e,i,j)*VVVO(b,c,e,k);
// 	//     X(a,b,c) -= t2(a,b,i,m)*OVOO(m,c,j,k);
// 	//     X(a,b,c) += t2(b,e,j,k)*VVVO(a,c,e,i);
// 	//     X(a,b,c) -= t2(b,a,j,m)*OVOO(m,c,k,i);
// 	// }

// 	// for (int b = 0; b < N; ++b) {
// 	//     X(c,a,b) += t2(a,e,i,k)*VVVO(c,b,e,j);
// 	//     X(c,a,b) -= t2(a,c,i,m)*OVOO(m,b,k,j);
// 	//     X(c,a,b) += t2(c,e,k,i)*VVVO(a,b,e,j);
// 	//     X(c,a,b) -= t2(c,a,k,m)*OVOO(m,a,i,j);
// 	// }

// 	// for (int a = 0; a < N; ++a) {
// 	//     X(c,b,a) += t2(c,e,k,j)*VVVO(b,a,e,i);
// 	//     X(c,b,a) -= t2(c,a,k,m)*OVOO(m,b,j,i);
// 	// }

//     }

//     void triples(int i, int j, int k) {
// 	tensor::tensor<4> t2, VVVO, OVOO;
// 	tensor::tensor<3> X;
// 	triples(i, j, k, t2, VVVO, OVOO, X);
//     }

// }

// #endif // CC_TRIPLES_HPP
