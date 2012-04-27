#include <rysq/core.hpp>
#include "transpose.hpp"
#include "externals/cxx/sugar.hpp"
#include "externals/cxx/utility/permute.hpp"

using namespace rysq;


void Quartet<Shell>::initialize() {
    PACK(shells_, (&a_, &b_, &c_, &d_));

    L_ =  detail::L::add(*this);
    size_ = detail::size::multiply(*this);
    K_ = detail::K::multiply(*this);

    nc_ = 1;
    hybrid_ = 0;
    foreach (const Shell* shell, this->shells_) {
	nc_ *= shell->nc;
	hybrid_ += shell->is_hybrid();
    }
    // std::cout << *this << std::endl;
}


// Quartet<Shell> Quartet<Shell>::transposed(int transpose_mask) const {
//     UNPACK((const Shell *A, *B, *C, *D), this->shells_);

//     int shuffle_mask = Transpose::shuffle_mask(transpose_mask);
//     cxx::utility::permute(A, B, C, D, shuffle_mask);
//     // std::cout << shuffle_mask << "\n";
//     return Quartet<Shell>(const_cast<Shell&>(*A),
// 			  const_cast<Shell&>(*B),
// 			  const_cast<Shell&>(*C),
// 			  const_cast<Shell&>(*D));
// }
