#include "config.hpp"

#include "cc/cc.hpp"
#include "cc/tensor.hpp"
#include "cc/detail/map.hpp"
#include "cc/detail/permute.hpp"

// #include "cc/triples/host.hpp"
// #define CC_TRIPLES_THREAD_T3 Host::T3
// #include "cc/triples/thread.ijk.hpp"
// #undef CC_TRIPLES_THREAD_T3

// #ifdef HAVE_CUDA
// #include "cc/triples/device.hpp"
// #endif

#include <utility>
#include <stdexcept>

#include <boost/typeof/typeof.hpp>
#include <boost/numeric/ublas/adaptor.hpp>
// #include <boost/numeric/ublas/io.hpp>

#include <boost/multi_array.hpp>
#include <boost/multi_array/multi_array_ref.hpp>
#include <boost/progress.hpp>

#include "foreach.hpp"
#include "array/hdf5.hpp"

#include "blas.hpp"
#include "omp.hpp"
#include "runtime.hpp"
#include "utility/timer.hpp"

#include "boost/utility/profiler.hpp" 

namespace cc {
namespace triples {
namespace detail {
    

    struct Cache {
	struct range {
	    explicit range(int i) : start_(i), stop_(i+1) {}
	    range(int start, int stop) : start_(start), stop_(stop) {}
	    int start() const { return start_; }
	    int stop() const { return stop_; }
	    int size() const { return stop_ - start_; }
	private:
	    int start_, stop_;
	};
	struct no_such_key : std::runtime_error {
	    no_such_key(const char *what) : std::runtime_error(what) {}
	};
	template<typename I, typename J, typename K, typename L>
	void load(const std::string &key, const Array &a,
		  const I &i, const J &j, const K &k, const L &l) {
	    range ri(i), rj(j), rk(k), rl(l);
	    size_t start[] = { ri.start(), rj.start(), rk.start(), rl.start() };
	    size_t stop[] = { ri.stop(), rj.stop(), rk.stop(), rl.stop() };
	    data_[key].resize(boost::extents
			      [rl.size()][rk.size()]
			      [rj.size()][ri.size()]);
	    a.get(data_[key].data(), start, stop);
	}
	Tensor<4>& operator[](const char *key) {
	    if (!data_.count(key)) throw no_such_key(key);
	    return data_.find(key)->second;
	}
    private:
	std::map<std::string, Tensor<4> > data_;
    };

    struct T3 {
	typedef Cache::range range;

	T3(size_t n, size_t N)
	    : t_(boost::extents[N][N][N][n][n][n]) {}

	static int min(int a, range r) {
	    return std::min(a, r.stop()-1);
	}

	template<class F>
	size_t apply(F &f, Cache &cache, range ra, range rb, range rc) {
	    size_t ops = 0;
	    omp::task task;
#pragma omp parallel reduction(+:ops)
	    for (int it = -1,
		     next = task++,
		     c = rc.start(); c < rc.stop(); ++c) {
		for (int b = rb.start(); b <= min(c,rb); ++b) {
		    for (int a = ra.start(); a <= min(b,ra); ++a) {
			++it;
			if (it != next) continue;
			next = task++;
			
			typename Tensor<3>::reference t
			    = t_[c-rc.start()][b-rb.start()][a-ra.start()];
			ops += f(t, cache, ra, rb, rc, a, b, c);
		    }
		}
	    }
	    return ops;
	}

	size_t evaluate(range ra, range rb, range rc, Cache &cache) {
	    typedef Tensor<3>::reference T;
	    size_t ops = 0;
	    omp::task task;
#pragma omp parallel reduction(+:ops)
	    for (int it = -1,
		     next = task++,
		     c = rc.start(); c < rc.stop(); ++c) {
		for (int b = rb.start(); b <= min(c,rb); ++b) {
		    for (int a = ra.start(); a <= min(b,ra); ++a) {
			++it;
			if (it != next) continue;
			next = task++;

			T t = t_[c-rc.start()][b-rb.start()][a-ra.start()];
			assign(t, 0);
			ops += evaluate(ra, rb, rc, a, b, c, cache, t);
		    }
		}
	    }
	    return ops;
	}

	Tensor<6>& tensor() { return t_; }
	const Tensor<6>& tensor() const { return t_; }

    private:
	Tensor<6> t_;

    private:

	static
	size_t evaluate(range ra, range rb, range rc,
			int A, int B, int C,
			Cache &cache, Tensor<3>::reference &t) {

	    int a = A - ra.start();
	    int b = B - rb.start();
	    int c = C - rc.start();

	    size_t ops = 0;

	    // t(i,j,k) = t(i,j,e,a) V(e,k,b,c) - t(i,m,a,b) V(j,k,m,c)
	    // t(i,k,j) = t(i,k,e,a) V(e,j,c,b) - t(i,m,a,c) V(k,j,m,b)
	    ops +=
		evaluate(cache["t(o,o,v,a)"][a], cache["V(v,o,v,c)"][c][B],
			 cache["t(o,o,v,a)"][a], cache["V(v,o,v,b)"][b][C],
			 cache["t(o,o,v,a)"][a][B], cache["V(o,o,o,c)"][c],
			 cache["t(o,o,v,a)"][a][C], cache["V(o,o,o,b)"][b],
			 t);

 	    // t(k,i,j) = t(k,i,e,c) V(e,j,a,b) - t(k,m,c,a) V(i,j,m,b)
	    // t(k,j,i) = t(k,j,e,c) V(e,i,b,a) - t(k,m,c,b) V(j,i,m,a)
	    ops +=
	    	evaluate(cache["t(o,o,v,c)"][c], cache["V(v,o,v,b)"][b][A],
	    		 cache["t(o,o,v,c)"][c], cache["V(v,o,v,a)"][a][B],
	    		 cache["t(o,o,v,c)"][c][A], cache["V(o,o,o,b)"][b],
	    		 cache["t(o,o,v,c)"][c][B], cache["V(o,o,o,a)"][a],
	    		 t);

	    // // t(j,k,i) = t(j,k,e,b) V(e,i,a,c) - t(j,m,b,c) V(k,i,m,a)
	    // // t(j,i,k) = t(j,i,e,b) V(e,k,c,a) - t(j,m,b,a) V(i,k,m,c)
	    ops +=
	    	evaluate(cache["t(o,o,v,b)"][b], cache["V(v,o,v,a)"][a][C],
	    		 cache["t(o,o,v,b)"][b], cache["V(v,o,v,c)"][c][A],
	    		 cache["t(o,o,v,b)"][b][C], cache["V(o,o,o,a)"][a],
	    		 cache["t(o,o,v,b)"][b][A], cache["V(o,o,o,c)"][c],
	    		 t);

	    return ops;
	}

	static
	size_t evaluate(Tensor<3>::reference t_ijae, Tensor<2>::reference V_ekbc,
			Tensor<3>::reference t_ikae, Tensor<2>::reference V_ejcb,
			Tensor<2>::reference t_imab, Tensor<3>::reference V_jkmc,
			Tensor<2>::reference t_imac, Tensor<3>::reference V_kjmb,
			Tensor<3>::reference t) {
	    using boost::numeric::ublas::trans;
	    size_t ops = 0;

	    ops += contract(1, as_matrix<2,1>(t_ijae), as_matrix<1,1>(V_ekbc),
	    		    as_matrix<2,1>(t));
	    ops += contract(-1, as_matrix<1,1>(t_imab), trans(as_matrix<2,1>(V_jkmc)),
	    		    as_matrix<1,2>(t));
	    permute<0,2,1>(t);

	    ops += contract(1, as_matrix<2,1>(t_ikae), as_matrix<1,1>(V_ejcb),
	    		    as_matrix<2,1>(t));
	    ops += contract(-1, as_matrix<1,1>(t_imac), trans(as_matrix<2,1>(V_kjmb)),
	    		    as_matrix<1,2>(t));
	    permute<1,0,2>(t);

	    return ops;
	}	
 
 	static void assign(Tensor<3>::reference a, double value) {
	    size_t size = 1;
	    for (size_t i = 0; i < 3; ++i) {
		size *= a.shape()[i];
	    }
	    std::fill(a.origin(), a.origin()+size, value);
	}

	template<class A, class B, class C>
	static size_t contract(double alpha, const A &a, const B &b, C c) {
	    assert(a.size1() == c.size1());
	    assert(b.size2() == c.size2());
	    assert(a.size2() == b.size1());
	    blas::set_num_threads(1);
	    blas::gemm(alpha, a, b, 1, c);
	    return 2*(a.size1()*a.size2()*b.size2());
	}

    };


    struct Energy  {
	const size_t no, nv;
	const Vector eh, ep;
    private:
	Matrix t1_, u_;
	Correction corr_;
	struct Partial {
	    Correction corr;
	    Vector u[3];
	    explicit Partial(size_t n) {
		foreach (Vector& v, u) {
		    v.resize(n);
		    v.clear();
		}
	    }
	};
	template<class T>
	struct V { T ab, ac, ba, bc, ca, cb; };

    public:
	Energy(const Array &t1, const Vector &eh, const Vector &ep) :
	    no(eh.size()), nv(ep.size()), eh(eh), ep(ep)
	{
	    t1_.resize(no, nv, false);
	    u_.resize(no, nv, false);
	    u_.clear();
	    size_t start[] = { 0, 0 };
	    size_t stop[] = { no, nv };
	    t1.get(t1_.data().begin(), start, stop);
	}

	Correction corr() const {
	    double ets = 0;
	    for (size_t i = 0; i < (nv*no); ++i) {
	    	ets += t1_.data()[i]*u_.data()[i];
	    }
	    Correction corr = corr_;
	    corr.ets += 2*ets;
	    return corr;
	}

	template<class T3>
	size_t operator()(const T3 &t3, Cache &cache,
			  Cache::range ra, Cache::range rb, Cache::range rc,
			  int a, int b, int c) {
	    double S = 1.0/(1 + (a == b || b == c));
	    S *= !(a == b && b == c); // zero diagonal
	    double D = ep(a) + ep(b) + ep(c);

	    Partial p = evaluate(no, S, D,
				 cache["V(o,o,v,b)"][b-rb.start()],
				 cache["V(o,o,v,c)"][c-rc.start()],
				 t3, eh, a, b, c);
#pragma omp critical
	    {
		using boost::numeric::ublas::column;
		column(u_, a) += p.u[0];
		column(u_, b) += p.u[1];
		column(u_, c) += p.u[2];
		corr_ += p.corr;
	    }
	    return 0;
	}

    private:
	template<class V2, class T3>
	static
	Partial evaluate(size_t no, double S, double Dp,
			 const V2 &Vb, const V2 &Vc,
			 const T3 &t3, const Vector &eh,
			 int a, int b, int c) {


#define t3(i,j,k) t3[k][j][i]
#define t(i,j,k) t ## i ## j ## k
#define V(i,j,a,b) (V ## b [a][j][i])

	    Partial p(no);
	    int n = no;
	    for (int k = 0; k < n; ++k) {
		for (int j = 0; j < n; ++j) {
		    double djk = eh(j) + eh(k);

		    struct { double bc, ac, ab; }
		    Vjk = { V(j,k,b,c), V(j,k,a,c), V(j,k,a,b) };
		    struct { double bc, ac, ab; }
		    Vkj = { V(k,j,b,c), V(k,j,a,c), V(k,j,a,b) };

		    for (int i = 0; i < n; ++i) {
			double D = S/(eh(i) + djk - Dp);

			double t(i,j,k) = t3(i,j,k);
			double t(i,k,j) = t3(i,k,j);
			double t(j,i,k) = t3(j,i,k);
			double t(j,k,i) = t3(j,k,i);
			double t(k,i,j) = t3(k,i,j);
			double t(k,j,i) = t3(k,j,i);

			double f = (8*t(i,j,k) -
				    4*(t(i,k,j)+t(k,j,i)+t(j,i,k)) +
				    2*(t(j,k,i)+t(k,i,j)));
			p.corr.etd += D*f*t(i,j,k);

			// ET[S]

			double Z[] = {
			    (2*t(i,j,k) + t(k,i,j)) - (2*t(j,i,k) + t(i,k,j)), // abc
			    (2*t(i,k,j) + t(k,j,i)) - (2*t(j,k,i) + t(i,j,k)), // acb
			    (2*t(j,i,k) + t(i,k,j)) - (2*t(i,j,k) + t(k,i,j)), // bac
			    (2*t(k,i,j) + t(j,k,i)) - (2*t(k,j,i) + t(j,i,k)), // bca
			    (2*t(j,k,i) + t(i,j,k)) - (2*t(i,k,j) + t(k,j,i)), // caa
			    (2*t(k,j,i) + t(j,i,k)) - (2*t(k,i,j) + t(j,k,i))  // cba
			};

			p.u[0](i) += D*(Z[0]*Vjk.bc);
			p.u[0](i) += D*(Z[1]*Vkj.bc);
			p.u[1](i) += D*(Z[2]*Vjk.ac);
			p.u[1](i) += D*(Z[3]*Vkj.ac);
			p.u[2](i) += D*(Z[4]*Vjk.ab);
			p.u[2](i) += D*(Z[5]*Vkj.ab);

		    }
		}
	    }
#undef t3
#undef t
#undef V
	    return p;
	}
    };

    struct Threads {

	explicit Threads(const runtime &rt) {}

	void operator()(size_t no, size_t nv,
			  const Map<const Array*> &t,
			  const Map<const Array*> &V,
			  const Map<Vector> &e) {

	    struct { size_t a, b, c; } last = { -1, -1, -1 };

	    const size_t block = 5;

	    std::cout << "Memory needed: "
	    	      << (block*3*(nv*nv*no + nv*no*no + no*no*no) +
			  block*3*(nv*no*no) +
	    		  block*block*block*(no*no*no))*8*1e-6
	    	      << " MB" << std::endl;

	    Cache cache;
	    T3 t3(no, block);

	    size_t nb = (nv + block - 1)/block;
	    boost::progress_display progress((nb*(nb+1)*(nb+2))/6);

	    utility::timer time;
	    boost::utility::timer::value_type io;

	    Energy energy(*t("i,a"), e("i"), e("a"));

	    size_t ops = 0;

	    using std::min;

	    for (size_t c = 0; c < nv; c += block) {
	    	for (size_t b = 0; b <= c; b += block) {
	    	    for (size_t a = 0; a <= b; a += block) {

	    		Cache::range v(0,nv), o(0,no);
	    		Cache::range ra(a, min(a + block, nv));
	    		Cache::range rb(b, min(b + block, nv));
	    		Cache::range rc(c, min(c + block, nv));

	    		boost::utility::timer timer;

	    		if (a != last.a) {
			    cache.load("V(v,o,v,a)", *V("e,k,b,c"), v, o, v, ra);
	    		    cache.load("t(o,o,v,a)", *t("i,j,a,b"), o, o, v, ra);
	    		    cache.load("V(o,o,o,a)", *V("j,k,i,a"), o, o, o, ra);
	    		}

	    		if (b != last.b) {
	    		    cache.load("t(o,o,v,b)", *t("i,j,a,b"), o, o, v, rb);
	    		    cache.load("V(o,o,o,b)", *V("j,k,i,a"), o, o, o, rb);
	    		    cache.load("V(v,o,v,b)", *V("e,k,b,c"), v, o, v, rb);
			    cache.load("V(o,o,v,b)", *V("i,j,a,b"), o, o, v, rb);
	    		}

	    		if (c != last.c) {
	    		    cache.load("t(o,o,v,c)", *t("i,j,a,b"), o, o, v, rc);
	    		    cache.load("V(o,o,o,c)", *V("j,k,i,a"), o, o, o, rc);
	    		    cache.load("V(v,o,v,c)", *V("e,k,b,c"), v, o, v, rc);
			    cache.load("V(o,o,v,c)", *V("i,j,a,b"), o, o, v, rc);
	    		}

	    		io += timer;

	    		ops += t3.evaluate(ra, rb, rc, cache);
			t3.apply(energy, cache, ra, rb, rc);

	    		last.a = a;
	    		last.b = b;
	    		last.c = c;

			++progress;
	    	    }
	    	}
	    }

	    Correction corr = energy.corr();

	    std::cout << "CCSD[T]: " << corr.etd << std::endl;
	    std::cout << "CCSD(T): " << corr.ets << std::endl;

	    std::cout << "I/O time: " << io << std::endl;
	    std::cout << "Total time: " << time << std::endl;
	    std::cout << "GFLOP/s: " << ops/1e9 << "/" << double(time)
	    	      << " = " << (ops/1e9)/time << std::endl;

	    return;
	}
    };

} // namespace triples
} // namespace detail
} // namespace cc


cc::Triples::Correction
cc::Triples::operator()(size_t no, size_t nv,
			Map<const Array*> t,
			Map<const Array*> V,
			const Map<Vector> &e) {

    //std::cout << e << std::endl;
    std::cout << "computing (t)" << std::endl;

    boost::utility::global_profiler().clear();

    triples::detail::Threads threads(runtime("CC"));
    threads(no, nv, t, V, e);


    Correction C; // = result.C;
    // triples::detail::Matrix t1(nv, no, 0);
    // const size_t r1[] = { 0, 0 };
    // const size_t r2[] = { nv, no };
    // t("a,i")->get(t1.data().begin(), r1, r2);
    // double ets = 0;
    // for (int i = 0; i < int(nv*no); ++i) {
    // 	ets += t1.data()[i]*result.u1.data()[i];
    // }
    // C.ets += 2*ets;

    std::cout << boost::utility::global_profiler() << std::endl;

    return C;

}
