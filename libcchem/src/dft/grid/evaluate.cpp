#include "dft/grid.hpp"
#include "blas.hpp"
#include "dft/wavefunction.hpp"
#include "dft/matrix.hpp"

#include <boost/numeric/ublas/vector.hpp>

namespace dft {

    template<class C, class Q, class O, class T>
    std::vector<double>
    gradient(const C &c, const Q &dq, const O &mo, T &tmp, double cutoff) {
	namespace blas = boost::numeric::bindings::blas;
	namespace ublas = boost::numeric::ublas;

	const int N = dq.size2();
	BOOST_AUTO(dmo, ublas::subrange(tmp, 0, c.size1(), 0, N)); 
	block::prod<16>(1.0, c, dq, 0.0, dmo, cutoff);
	//blas::gemm(1.0, c, dq, 0.0, dmo);

	std::vector<double> r(N);
#pragma omp parallel for schedule(dynamic,16)
	for (int i = 0; i < N; ++i) {
	    r[i] = blas::dot(ublas::column(mo,i), ublas::column(dmo,i));
	    r[i] *= 2;
	}

	return r;
    }

    void Grid::evaluate(const Wavefunction &W,
			const_array_ref<double,3> dr,
			const_array_ref<double> w) {

	namespace ublas = boost::numeric::ublas;
	namespace blas = boost::numeric::bindings::blas;
	BOOST_AUTO(const &basis, W.basis());

	BOOST_AUTO(const &Ct, W.C(W.occupied()));
	BOOST_AUTO(C, ublas::trans(Ct));

	const int N = dr.size();
	BOOST_AUTO(mo, ublas::project(vectors.mo(),
				      ublas::range(0, C.size1()),
				      ublas::range(0, N)));

#pragma omp parallel
	{
	    ublas::vector<double> ao(basis.size());
#pragma omp for schedule(dynamic,16)
	    for (int i = 0; i < N; ++i) {
		dft::wavefunction::vector<0> f0(&ao[0]);
		dft::wavefunction::evaluate<0>(basis, Point(dr[i]), f0);
		vectors.ao(i) = ao;
	    }
	}

	block::prod<16>(1.0, C, vectors.ao(), 0.0, mo, p_.ccut);
	//blas::gemm(1, C, vectors.ao(), 0, mo);

#pragma omp parallel for schedule(dynamic,16) ordered
	for (int i = 0; i < N; ++i) {
	    double rho[2] = { 0 };
	    for (int s = 0; s < 1; ++s) {
		BOOST_AUTO(const &v, vectors.mo(i));
		rho[s] = blas::dot(v,v);
		if (rho[s] < 0) rho[s] = 0;
	    }
	    if ((rho[0] + rho[1]) < p_.rcut) continue;
#pragma omp ordered
#pragma omp critical
	    density.push_back(i, w[i], rho, 1);
	}

	for (size_t i = 0; i < density.size(); ++i) {
	    size_t j = density.index[i];
	    if (i != j) {
		vectors.ao(i) = vectors.ao(j);
		for (int s = 0; s < 1; ++s) {
		    vectors.mo(i) = vectors.mo(j);
		}
	    }
	}
    }


    void Grid::evaluate(const Wavefunction &W,
			const_array_ref<double,3> dr) {

	namespace ublas = boost::numeric::ublas;
	namespace blas = boost::numeric::bindings::blas;
	BOOST_AUTO(const &basis, W.basis());

	BOOST_AUTO(const &Ct, W.C(W.occupied()));
	BOOST_AUTO(C, ublas::trans(Ct));

	const int N = density.size();

#pragma omp parallel
	{
	    ublas::vector<double> dx(basis.size());
	    ublas::vector<double> dy(basis.size());
	    ublas::vector<double> dz(basis.size());
#pragma omp for schedule(dynamic,16)
	    for (int i = 0; i < N; ++i) {
		dft::wavefunction::vector<1> f(&dx[0], &dy[0], &dz[0]);
		dft::wavefunction::evaluate<1>(basis, Point(dr[i]), f);
		vectors.dx(i) = dx;
		vectors.dy(i) = dy;
		vectors.dz(i) = dz;
	    }
	}

	BOOST_AUTO(&tmp, vectors.tmp());
	density.dx = dft::gradient(C, vectors.dx(), vectors.mo(), tmp, p_.ccut);
	density.dy = dft::gradient(C, vectors.dy(), vectors.mo(), tmp, p_.ccut);
	density.dz = dft::gradient(C, vectors.dz(), vectors.mo(), tmp, p_.ccut);

	for (int i = 0; i < N; ++i) {
	    density.aa.push_back(pow(density.dx[i],2) +
				 pow(density.dy[i],2) +
				 pow(density.dz[i],2));
	}
    }

} // namespace dft
