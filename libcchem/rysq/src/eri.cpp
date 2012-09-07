#include "rysq-eri.hpp"

#include "kernel/new.hpp"
#include "eri-transform.hpp"
#include "foreach.hpp"

#include <boost/typeof/typeof.hpp>

using namespace rysq;

typedef kernel::Eri<> kernel_type;

Eri::Eri(const Quartet<Shell> &quartet) {
    bool t[] = { quartet[0].size == 1 && quartet[1].size > 1,
		 quartet[2].size == 1 && quartet[3].size > 1,
		 quartet.bra().size() == 1 && quartet.ket().size() > 1 };
    transpose_ = Transpose(t[0], t[1], t[2]);
    kernel_ = kernel::new_<eri::Transform>(transpose_(quartet));
    size_ = quartet.size();
    // std::cout << transpose_(quartet) << std::endl;
}

Eri::~Eri() { delete static_cast<kernel_type*>(kernel_); }

void Eri::operator()(const Quartet<Center> &centers,
		     double *Q, const Parameters &parameters) {
    eri::Transform<>::Data data(Q);
    Parameters p = parameters;
    (*static_cast<kernel_type*>(kernel_))(transpose_(centers), data, p);
}

void Eri::operator()(const std::vector<Center> &centers,
		     const std::vector<Int4> &quartets,
		     double *Q, const Parameters &parameters) {
    // kernel_type &kernel = (*static_cast<kernel_type*>(kernel_));

    foreach (const Int4 &index, quartets) {
	(*this)(centers, index, Q, parameters);
	Q += size_;
    }
}

void Eri::operator()(const std::vector<Center> &centers, size_t size,
		     const Int4 quartets[], double* const Q[],
		     const Parameters &parameters) {
    for (size_t i = 0; i < size; ++i) {
	(*this)(centers, quartets[i], *(Q++), parameters);
    }
}
