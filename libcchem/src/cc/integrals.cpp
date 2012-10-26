// These headers are first - otherwise bizarre errors ith gcc < 4.3 and IMPI
#include "runtime.hpp"
#include "omp.hpp"
#include "blas.hpp"

#include "cc/cc.hpp"
#include "cc/tensor.hpp"
#include "cc/utility.hpp"

#include "thread.hpp"
#include "core/wavefunction.hpp"
#include "basis/basis.hpp"
#include "integrals/eri.hpp"
#include "foreach.hpp"
#include "utility/timer.hpp"
#include "utility/progress.hpp"


#include <boost/bind.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/adaptors.hpp>

#include "boost/utility/profiler.hpp"

namespace cchem {
namespace cc {
namespace detail {

    template<class R>
    void copy(const double *shell, boost::multi_array_ref<double,4> &a,
	      const R &p, const R &q, const R &r, const R &s) {
	for (size_t l = s.start(); l < s.stop(); ++l) {
	    for (size_t k = r.start(); k < r.stop(); ++k) {
		for (size_t j = q.start(); j < q.stop(); ++j) {
		    std::copy(shell, shell+p.size(), &a[j][k][l][p.start()]);
		    shell += p.size();
		}
	    }
	}
    }


    template<class C, class T>
    void transform12(const Basis &basis,
		     const Basis::Shell &R, const Basis::Block &S,
		     integrals::Eri eri, const C &C1, const C &C2,
		     T &t2, omp::lock &lock) {

	namespace ublas = boost::numeric::ublas;
	using tensor::as_matrix;
	typedef tensor_reference<4> T4;
	typedef typename T4::extent_range range;
	using boost::extents;

	const size_t MAX = std::min<size_t>(basis.max_block(), 128);
	BOOST_AUTO(const &blocks, basis.blocks(MAX));

	range rs = S.range();
	range rr = R.range();

	Thread thread;
	thread.data[0] = thread.malloc(MAX*(MAX*rr.size()*rs.size()));
	thread.data[1] = thread.malloc(C1.size2()*(MAX*rr.size()*rs.size()));

#pragma omp for schedule(dynamic,1)
	for (int q = 0; q < blocks.size(); ++q) {
	    BOOST_PROFILE_LINE;
	    utility::timer timer;

	    const Basis::Block &Q = blocks.at(q);

	    range rq = Q.range();
	    T4 t1(thread.data[1], extents[rq.size()][rr.size()][rs.size()][C1.size2()]);
	    std::fill_n(t1.data(), t1.num_elements(), 0);

	    foreach (const Basis::Block &P, blocks) {
		range rp = P.range();

		T4 G(thread.data[0], extents[rq][rr][rs][rp]);
		std::fill_n(G.data(), G.num_elements(), 0);

		eri(P, Q, R, S);
		size_t size = (P.shell().size()*Q.shell().size()*
			       R.size()*S.shell().size());
		for (size_t i = 0; i < eri.quartets().size(); ++i) {
		    BOOST_AUTO(const &q, eri.quartets()[i]);
		    detail::copy(eri.data() + i*size, G,
				 basis[q[0]], basis[q[1]],
				 basis[q[2]], basis[q[3]]);
		}
		ublas::range r1(rp.start(), rp.finish());
		ublas::range r2(0, C1.size2());
		BOOST_AUTO(t, (as_matrix<1,3>(t1)));
		blas::gemm(1,
			   project(trans(C1), r2, r1),
			   (as_matrix<1,3>(G)),
			   1, t);
	    }

	    // 2nd transformation
	    {
		BOOST_PROFILE_LINE;
		BOOST_AUTO(t, (as_matrix<3,1>(t2)));
		ublas::range r1(rq.start(), rq.finish());
		ublas::range r2(0, C2.size2());
		lock.set();
		blas::gemm(1,
			   tensor::as_matrix<3,1>(t1),
			   project(C2, r1, r2),
			   1, t);
		lock.unset();
	    }
	} // for q

    }


    // (O,N,s)
    template<class R, class C, class T2>
    void transform3(int q, R rs, const C &c, const T2 &t2,
		    Array &v, double *tmp) {
	BOOST_PROFILE_LINE;
	tensor_reference<3>
	    t3(tmp, boost::extents[rs.size()][c.size2()][t2.shape()[2]]);
	for (int i = 0; i < rs.size(); ++i) {
	    using tensor::as_matrix;
	    BOOST_AUTO(t, (as_matrix<1,1>(t3[i])));
	    blas::gemm(1, (as_matrix<1,1>(t2[i])), c, 0, t);
	}
	size_t start[4] = { 0, 0, q, rs.start() };
	size_t finish[4] = { t3.shape()[2], t3.shape()[1],
			     q+1, rs.start()+rs.size() };
	BOOST_PROFILE_LINE;
	v.put(t3.data(), start, finish);
    }


    template<class U, class C, class T>
    void transform(const U &u, const C &c, T &t) {
	// blas::gemm(1, u, c, 0, t);
	const int B = 128;
	using boost::numeric::ublas::range;
	int m = u.size1();
#pragma omp parallel for schedule(dynamic,1)
	for (int i = 0; i < m; i += B) {
	    range K(0, u.size2());
	    range M(i, std::min(i+B, m));
	    range N(0, c.size2());
	    BOOST_AUTO(ti, project(t, M, N));
	    blas::gemm(1, project(u, M, K), c, 0, ti);
	}
    }


    template<class T, class C>
    void transform4(int r, const C &c, const T &t, Array &V, double *buffer) {
	BOOST_PROFILE_LINE;
	{
	    using tensor::as_matrix;
	    using boost::numeric::ublas::make_matrix;
	    using boost::numeric::ublas::column_major;
	    int m = t.shape()[3]*t.shape()[2];
	    BOOST_AUTO(t4, make_matrix<column_major>(m, c.size2(), buffer));
	    transform(as_matrix<2,2>(t), c, t4);
	}
	size_t start[] = { 0, 0, r, 0 };
	size_t finish[] = { t.shape()[3], t.shape()[2], r+1, c.size2() };
	V.put(buffer, start, finish);
    }


}
}
} // namespace cchem


void cchem::cc::integrals(const Parallel &pe, const Wavefunction &wf,
			  const Map<Array*> &v,
			  const integrals::Screening &screening) {

    BOOST_PROFILE_REGISTER_THREAD;

    struct {
	utility::timer total;
	utility::timer::value_type eri;
	utility::timer::value_type trans[4];
    } time;

    namespace ublas = boost::numeric::ublas;
    using ublas::make_vector;
    using ublas::make_matrix;
    using ublas::column_major;
    using tensor::as_matrix;

    const Basis &basis = wf.basis();
    size_t N = basis.size();
    size_t no = wf.active().size();
    size_t nv = wf.virtuals().size();
    Matrix Co = wf.C(wf.active()); // C(p,i)
    Matrix Cv = wf.C(wf.virtuals()); // C(p,a)

    const size_t MAXS = basis.max().size();
    const size_t MAXR = basis.max().size();
    const size_t MAXP = std::min<size_t>(basis.max_block(), 128);
    const size_t MAXQ = MAXP;

    typedef tensor_reference<3> T3;
    typedef tensor_reference<4> T4;

    const bool IABC = v.has("iabc");

    typedef Symbol< tensor_reference<4> >::range range;
    range ro(0,no), rv(0,nv), ra(0,N);
    Matrix C(N,no+nv);
    subrange(C, 0, N, 0, nv) = Cv;
    subrange(C, 0, N, nv, nv+no) = Co;

    struct {
	std::vector<Basis::Block> P, Q, R, S;
    } blocks;

    // block shells into segments
    blocks.P = basis.blocks(MAXP);
    blocks.Q = basis.blocks(MAXQ);
    blocks.S = basis.blocks(MAXS);

    int max_threads = omp_get_max_threads();
    int nested = 1;
    {
	int threads = pe.size()*max_threads;
	int shells = basis.shells().size();
	nested = (threads+shells-1)/shells;
	nested = std::min(max_threads, nested);
    }

    {
	BOOST_AUTO(cout, pe.cout());
	cout << "Integral transformation: ";
	cout << ((IABC) ? "(iabc)" : "(ija*), (iajb)") << std::endl;
	cout << "Threads: " << max_threads << "/" << nested << std::endl;
    }

    utility::Progress progress;
    if (pe.rank() == 0) progress.reset(basis.size());

    foreach (const Basis::Block &S, blocks.S) {

	BOOST_PROFILE_LINE;

	pe.barrier();

	typedef T4::extent_range range;    
	range rs = S.range();

	omp_set_num_threads(max_threads/nested);

#pragma omp parallel
	{
	    BOOST_PROFILE_LINE;

	    Thread::Task<Parallel::Task&> task(pe.task());
	    task.reset();

	    Thread thread;
	    thread.buffer[2].resize(no*C.size2()*MAXR*MAXS);

	    while (task++ < basis.shells().size()) {
		BOOST_PROFILE_LINE;
		const Basis::Shell &R = basis[task];

		range rr = R.range();
		using boost::extents;
		T4 t2(thread.buffer[2], extents[C.size2()][rr.size()][rs.size()][no]);
		std::fill_n(t2.data(), t2.num_elements(), 0);

		omp_set_nested(nested > 1);
		omp_set_num_threads(nested);
		omp::lock lock;
#pragma omp parallel
		{
		    integrals::Eri eri(basis, &screening);
		    detail::transform12(basis, R, S, eri, Co, C, t2, lock);
		} // parallel

		{
		    BOOST_PROFILE_LINE;
		    size_t start[] = { 0, 0, rr.start(), 0 };
		    size_t finish[4];
		    for (int i = 0; i < 4; ++i) {
			finish[i] = start[i] + t2.shape()[3-i];
		    }
		    //std::cout <<  ns << " " << rr.start()+nr << " " << n2 << std::endl;
		    v["iqrs"]->put(t2.data(), start, finish);
		}

	    }
	} // parallel

	// restore threads
	omp_set_num_threads(max_threads);

	{
	    BOOST_PROFILE_LINE;
	    v["iqrs"]->flush();
	    pe.barrier();
	}

#pragma omp parallel
	{
	    BOOST_PROFILE_LINE;
	    
	    Thread::Task<Parallel::Task&> task(pe.task());
	    Thread thread;
	    BOOST_AUTO(&buffer, thread.buffer);
	    buffer[0].resize(no*rs.size()*N);
	    buffer[1].resize(no*N*rs.size());

	    if (!IABC) {
		task.reset();
		while (task++ < nv+no) {
		    int u = task;
		    T3 t2(buffer[0], boost::extents[rs.size()][N][no]);
		    {
			BOOST_PROFILE_LINE;
			T3 t(buffer[1], boost::extents[N][rs.size()][no]);
			size_t start[] = { 0, 0, 0, u };
			size_t finish[] = { no, rs.size(), N, u+1 };
			v["iqrs"]->get(t.data(), start, finish);
			tensor::copy<0,2,1>(t, t2);
		    }
		    using detail::transform3;
		    // (O,VO,q,i) -> (O,O,VO,s)
		    transform3(u, rs, Co, t2, *v["ijab"], buffer[1]);
		    // (O,s,q,u) -> (O,V,u,s)
		    if (u < nv) continue;
		    transform3(u-nv, rs, Cv, t2, *v["iajb"], buffer[1]);
		}
	    }

	    // (O,V,r,s) -> (O,rs,V)
	    if (IABC) {
		BOOST_PROFILE_LINE;
		task.reset();
		while (task++ < nv) {
		    size_t u = task;
		    T3 t2(buffer[0], boost::extents[N][rs.size()][no]);
		    {
			size_t start[] = { 0, 0, 0, u };
			size_t finish[] = { no, rs.size(), N, u+1 };
			v["iqrs"]->get(t2.data(), start, finish);
		    }
		    T3 t3(buffer[0], boost::extents[rs.size()][nv][no]);
		    {
			T3 t_(buffer[1], boost::extents[nv][rs.size()][no]);
			BOOST_AUTO(t, (as_matrix<2,1>(t_)));
			blas::gemm(1, (as_matrix<2,1>(t2)), Cv, 0, t);
			tensor::copy<0,2,1>(t_, t3);
		    }
		    size_t start[] = { 0, 0, rs.start(), u };
		    size_t finish[] = { no, nv, rs.finish(), u+1 };
		    v["iabc"]->put(t3.data(), start, finish);
		}
	    }

	} // parallel

	progress += rs.size();

    } // foreach S

    {
	BOOST_PROFILE_LINE;
	if (!IABC) {
	    v["iajb"]->flush();
	    v["ijab"]->flush();
	}
	if (IABC) {
	    v["iabc"]->flush();
	}
	pe.barrier();
    }

    BOOST_PROFILE_DUMP(std::cout);

    if (!IABC) {
	BOOST_PROFILE_LINE;

	Thread::Task<Parallel::Task&> task(pe.task());

	using detail::transform4;

	Symbol< tensor_reference<4> > S;

	Buffer<double> buffer;
	buffer[0].resize(no*nv*N);
	buffer[1].resize(no*nv*nv);

	task.reset();
	while (task++ < nv) {
	    int i = task;
	    S.load("ijab", buffer[0], ro, ro, i, ra, *v["ijab"]);
	    transform4(i, Cv, S["ijab"], *v["ijab"], buffer[1]);
	}

	task.reset();
	while (task++ < no) {
	    int i = task;
	    S.load("ijkf", buffer[0], ro, ro, i+nv, ra, *v["ijab"]);
	    transform4(i, Cv, S["ijkf"], *v["ijka"], buffer[1]);
	    transform4(i, Co, S["ijkf"], *v["ijkl"], buffer[1]);
	}

	task.reset();
	while (task++ < no) {
	    int i = task;
	    S.load("iajb", buffer[0], ro, rv, i, ra, *v["iajb"]);
	    transform4(i, Cv, S["iajb"], *v["iajb"], buffer[1]);
	}

	BOOST_PROFILE_LINE;
	v["ijka"]->flush();
	v["ijab"]->flush();
	v["iajb"]->flush();
	v["ijab"]->flush();
	pe.barrier();
    }

    // 
    if (IABC) {
	BOOST_PROFILE_LINE;
	using detail::transform4;
	utility::timer timer;
	Symbol< tensor_reference<4> > S;

	Buffer<double> buffer;
	buffer[0].resize(no*nv*N);
	buffer[1].resize(no*nv*nv);

	Thread::Task<Parallel::Task&> task(pe.task());
	task.reset();

	while (task++ < nv) {
	    int c = task;
	    S.load("t3", buffer[0], ro, rv, ra, c, *v["iabc"]);
	    // 3-index
	    S.set("t4", buffer[1], ro, rv, rv, c);
	    {
		BOOST_AUTO(t4, (tensor::as_matrix<2,1>(S["t4"][0]))); 
		blas::gemm(1, tensor::as_matrix<2,1>(S["t3"][0]), Cv, 0, t4);
	    }
	    S.store("t4", ro, rv, rv, c, *v["iabc"]);
	}
	time.trans[2] = timer;

	BOOST_PROFILE_LINE;
	v["iabc"]->flush();
	pe.barrier();
    }	

    BOOST_PROFILE_DUMP(std::cout);
	
}


