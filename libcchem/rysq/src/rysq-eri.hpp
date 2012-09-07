#ifndef RYSQ_ERI_HPP
#define RYSQ_ERI_HPP

#include "rysq-core.hpp"

#include <vector>
#include <memory>

namespace rysq {

    class Eri {
    public:
	class Impl;
	size_t size;

	struct Parameters {
	    double cutoff, scale;
	    Parameters(double cutoff = 1.0e-10, double scale = 1)
		: cutoff(cutoff), scale(scale) {}
	};

	Eri(const Quartet<Shell> &shells);
	~Eri();
	void operator()(const Quartet<Center> &centers, double *I,
			const Parameters &parameters = Parameters());

	void operator()(const std::vector<Center> &centers,
			const Int4 &quartet, double *I,
			const Parameters &parameters = Parameters()) {
	    Quartet<Center> r(centers.at(quartet[0]), centers.at(quartet[1]),
			      centers.at(quartet[2]), centers.at(quartet[3]));
	    (*this)(r, I, parameters);
	}

	void operator()(const std::vector<Center> &centers,
			const std::vector<Int4> &quartets,
			double *Eri,
			const Parameters &parameters = Parameters());

	void operator()(const std::vector<Center> &centers, size_t size,
			const Int4 quartets[], double* const Eri[],
			const Parameters &parameters = Parameters());
    private:
	Transpose transpose_;
	void *kernel_;
	size_t size_;
    };

}


#endif /* RYSQ_ERI_HPP */
