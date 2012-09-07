#ifndef RYSQ_TRANSPOSE_HPP
#define RYSQ_TRANSPOSE_HPP

#include "rysq-core.hpp"
// #include "externals/cxx/utility/permute.hpp"

namespace rysq {


// inline int transpose_permutation(int value) {
//     static const int index[2][2][2] = {{{ SHUFFLE(0,1,2,3), SHUFFLE(1,0,2,3) },
// 					{ SHUFFLE(0,1,3,2), SHUFFLE(1,0,3,2) }},
// 				       {{ SHUFFLE(2,3,0,1), SHUFFLE(2,3,1,0) },
// 					{ SHUFFLE(3,2,0,1), SHUFFLE(3,2,1,0) }}};
//     bool bra = value & Transpose::BRA;
//     bool ket = value & Transpose::KET;
//     bool braket = value & Transpose::BRAKET;
//     // std::cout << bra << ket << braket  << " " << index[bra][ket][braket] << std::endl;
//     return index[braket][ket][bra];
// }

template<class Q> Q transpose(Q q, int value) {
    if (value & Transpose::BRA) std::swap(q[0], q[1]);
    if (value & Transpose::KET) std::swap(q[2], q[3]);
    if (value & Transpose::BRAKET) {
	std::swap(q[0], q[2]);
	std::swap(q[1], q[3]);
    }
    return q;
}

inline Quartet<Shell> transpose(const Quartet<Shell> &quartet, int value) {
    boost::array<Shell*,4> p = transpose(quartet.data(), value);
    // std::cout << p << quartet.data() << std::endl;
    return Quartet<Shell>(const_cast<Shell&>(*p[0]),
			  const_cast<Shell&>(*p[1]),
			  const_cast<Shell&>(*p[2]),
			  const_cast<Shell&>(*p[3]));
}

template<class Q> Q Transpose::operator()(const Q &quartet) {
    return transpose(quartet, this->value);
}


namespace  mpl {

    template<bool bra, bool ket, bool braket>
    struct transpose {
	static const int value =
	    ((int(bra) << 0) | (int(ket) << 1) | (int(braket) << 2));
    };

}

template<int ni, int nj, int nk, int nl>
void transpose(int value, double *Q) {

#define TRANSPOSE_(I,J,K,L) {						\
	double T[n ## L][n ## K][n ## J][n ## I] __attribute__ ((aligned(16))); \
	for (int l = 0, ijkl = 0; l < nl; ++l) {			\
	    for (int k = 0; k < nk; ++k) {				\
		for (int j = 0; j < nj; ++j) {				\
		    for (int i = 0; i < ni; ++i, ++ijkl) {		\
			T[L][K][J][I] = Q[ijkl];			\
		    }							\
		}							\
	    }								\
	}								\
	static const int N = ni*nj*nk*nl;				\
	for (int i = 0; i < N; ++i)					\
	    Q[i] = (&(T[0][0][0][0]))[i];				\
    }

    static const int BRA = Transpose::BRA;
    static const int KET = Transpose::KET;
    static const int BRAKET = Transpose::BRAKET;

    if (value & BRAKET) {
	if (value == (BRA | KET | BRAKET)) { TRANSPOSE_(l,k,j,i); }
	else if (value == (BRA | BRAKET))  { TRANSPOSE_(l,k,i,j); }
	else if (value == (KET | BRAKET))  { TRANSPOSE_(k,l,j,i); }
	else                                   { TRANSPOSE_(k,l,i,j); }
    }
    else {
	if (value == (BRA | KET)) { TRANSPOSE_(j,i,l,k); }
	else if (value == (BRA))  { TRANSPOSE_(j,i,k,l); }
	else if (value == (KET))  { TRANSPOSE_(i,j,l,k); }
    }	

#undef TRANSPOSE_

}

} // namespace rysq

#endif // RYSQ_TRANSPOSE_HPP
