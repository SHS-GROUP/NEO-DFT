#include <boost/thread/tss.hpp>
#include <boost/typeof/typeof.hpp>

#include "rysq-cuda.hpp"
#include "cuda/fock.hpp"
#include "transpose.hpp"
#include "foreach.hpp"

#include "boost/utility/profiler.hpp"

rysq::cuda::Fock::Fock(Shell &a, Shell &b, Shell &c, Shell &d,
		       const Context &context) {

    bool t0 = (!a.L() && b.L());
    bool t1 = (!c.L() && d.L());
    bool t01 = (!(a.L() + b.L()) && (c.L() + d.L()));
    Transpose transpose(t0, t1, t01);
    detail::Shell::Quartet quartet(a, b, c, d);

    try {
	this->kernel_.reset(new cuda::detail::Fock(quartet, context, transpose));
	for (int i = 0; i < 4; ++i)
	    block_[i] = quartet[i].size();
    }
    catch (...) {
	this->kernel_.reset();
    }
}

rysq::cuda::Fock::~Fock() {}

void rysq::cuda::Fock::operator()(const Centers &centers,
				  const Quartets &quartets,
				  density_matrix_set D, fock_matrix_set F,
				  Mutex &mutex,
				  const Parameters &parameters,
				  Stream s) {
    detail::matrix_set<double> S(D, F);
    (*kernel_)(detail::Centers(centers),
	       detail::Quartets(quartets),
	       S, detail::Mutex(mutex), parameters,
	       static_cast<cudaStream_t>(s.data()));
    // {
    // 	BOOST_PROFILE_LINE;
    // 	s.synchronize();
    // }
}
