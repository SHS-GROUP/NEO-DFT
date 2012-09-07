#include "integrals/screening.hpp"
#include "integrals/rysq.hpp"
#include "integrals/eri.hpp"
#include "basis/basis.hpp"

namespace integrals {

    void initialize() {
	integrals::rysq::initialize();
    }

    Screening::Screening(const Basis &basis, double value) {
	int N = basis.shells().size();

	K_.resize(N,N);
	K_.clear();
	value_ = value;

	// std::cout << value_ << std::endl;;

	Screening S(N);
	integrals::Eri eri(basis, &S);

	for (int j = 0; j < N; ++j) {
	    for (int i = 0; i <= j; ++i) {
		const Basis::Shell &P = basis.shells().at(i);
		const Basis::Shell &Q = basis.shells().at(j);
		eri(P,Q,P,Q);
		double max = 0;
		size_t n = P.size()*Q.size();
		if (!eri.quartets().empty()) {
		    for (size_t k = 0; k < n*n; ++k) {
			max = std::max(max, fabs(eri.data()[k]));
		    }
		}
		K_(i,j) = sqrt(max);
		K_(j,i) = sqrt(max);
	    }
	}
    }

} // namespace integrals
