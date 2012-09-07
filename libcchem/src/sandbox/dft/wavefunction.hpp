#ifndef DFT_WAVEFUNCTION_HPP
#define DFT_WAVEFUNCTION_HPP

#include "foreach.hpp"

#include "basis/config.hpp"
#include "basis/basis.hpp"
#include "fusion/shell.hpp"

// #include <vector>

#include <boost/preprocessor/seq/for_each.hpp>

// #include <boost/array.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/static_assert.hpp>

#include <boost/fusion/include/algorithm.hpp>

namespace dft {
namespace wavefunction {

    template<int N, class I = double*>
    struct vector;

    template<class I>
    struct vector<0,I> {
	typedef typename std::iterator_traits<I>::reference reference;
	I data;
	vector(I data) : data(data) {}
	reference operator*() { return *data; }
	void operator++() { ++data; }
	vector operator+(int dist) const {
	    return vector(data + dist);
	}
	void clear_range(int begin, int end) {
	    for (int i = begin; i < end; ++i) {
		*(data+i) = 0;
	    }
	}
    };

    template<class I>
    struct vector<1,I> {
	I dx, dy, dz;
	vector(I dx, I dy, I dz) : dx(dx), dy(dy), dz(dz) {}
	void operator++() { ++dx; ++dy; ++dz; }
	vector operator+(int dist) const {
	    return vector(dx + dist, dy + dist, dz + dist);
	}
	void clear_range(int begin, int end) {
	    for (int i = begin; i < end; ++i) {
		*(dx+i) = *(dy+i) = *(dz+i) = 0;
	    }
	}
    };


    struct shell {

	template<int N>
	struct factors {
	    static const int size = (N) ? N : 1;
	    double data[size][3];

	    template<class A>
	    factors(const A &r) {
		for (int q = 0; q < 3; ++q) {
		    data[0][q] = r[q];
		    // std::cout << data[0][q] << std::endl;
		    for (int j = 1; j < size; ++j) {
			data[j][q] = data[0][q]*data[j-1][q];
		    }
		}
	    }

	    static double pow2(const double &x) { return x*x; }

	    double r2() const {
		return pow2(get<0,1>()) + pow2(get<1,1>()) + pow2(get<2,1>());
	    }

	    template<int l, int m, int n>
	    typename boost::disable_if_c<
	    	(l < 0 || m < 0 || n < 0), double>::type
	    xyz() const {
		return get<0,l>()*get<1,m>()*get<2,n>();
	    }

	    template<int l, int m, int n>
	    typename boost::enable_if_c<
	    	(l < 0 || m < 0 || n < 0), int>::type
	    xyz() const { return 0; }

	    template<int Q, int L>
	    typename boost::enable_if_c<(L > 0), const double&>::type
	    get() const {
		// BOOST_STATIC_ASSERT((0 < L  && L <= size));
		return data[L-1][Q];
	    }
            template<int Q, int L>
	    typename boost::disable_if_c<(L > 0), double>::type
	    get() const { return (L == 0); }
	};

	template<int L, int N>
	struct primitive;

	template<int L>
	struct primitive<L,0> {
	    double c;
	    const factors<L> &q;
	    mutable vector<0> v;
	    primitive(double C, double alpha, const factors<L> &q, vector<0> v)
		: c(C*exp(-alpha*q.r2())), q(q), v(v) {}
	    bool test(double value) const { return (c > value); }
	    template<class F>
	    void operator()(const F&) const {
		*v += c*q.template xyz<F::l, F::m, F::n>();
		++v;
	    }
	};

	template<int L>
	struct primitive<L,1> {
	    double c, a;
	    const factors<L+1> &q;
	    mutable vector<1> v;
	    primitive(double C, double alpha, const factors<L+1> &q, vector<1> v)
		: c(C*exp(-alpha*q.r2())), a(2*alpha), q(q), v(v) {}
	    bool test(double value) const { return (c > value); }
	    template<class F>
	    void operator()(const F &f) const {
		*v.dx += c*(q.template xyz<F::l-1, F::m, F::n>()*F::l
			    - a*(q.template xyz<F::l+1, F::m, F::n>()));
		*v.dy += c*(q.template xyz<F::l, F::m-1, F::n>()*F::m
			    - a*(q.template xyz<F::l, F::m+1, F::n>()));
		*v.dz += c*(q.template xyz<F::l, F::m, F::n-1>()*F::n
			    - a*(q.template xyz<F::l, F::m, F::n+1>()));
		++v;
	    }
	};

	template<basis::shell::type type, int N, class R>
	static void
	evaluate(const basis::Block &block, const R &r, vector<N> v) {
	    static const int L = cchem::fusion::shell<type>::L;
	    for (BOOST_AUTO(shell, block.begin()); shell < block.end(); ++shell) {
		factors<L+N> q(r(*shell));
		evaluate<type>(*shell, q, v + shell->start());
	    }
	}

    private:

	template<basis::shell::type type, int M, int N>
	static typename boost::enable_if_c<(type > -1)>::type
	evaluate(const basis::Shell::Data &shell,
		 const factors<M> &q, vector<N> v) {
	    typename cchem::fusion::shell<type>::functions f;
	    for (int k = 0; k < shell.K(); ++k) {
		primitive<type,N> p(shell.C(0)[k], shell.a()[k], q, v);
		//if (!p.test(1e-10)) continue;
		boost::fusion::for_each(f, p);
	    }
	}

	template<basis::shell::type type, int M, int N>
	static typename boost::enable_if_c<(type == -1)>::type
	evaluate(const basis::Shell::Data &shell,
		 const factors<M> &q, vector<N> v) {
	    static const int L = 1;
	    typename cchem::fusion::shell<basis::shell::S>::functions S;
	    typename cchem::fusion::shell<basis::shell::P>::functions P;
	    for (int k = 0; k < shell.K(); ++k) {
		primitive<L,N> s(shell.C(0)[k], shell.a()[k], q, v);
		boost::fusion::for_each(S, s);
		primitive<L,N> p(shell.C(1)[k], shell.a()[k], q, v+1);
		//if (!p.test(1e-10)) continue;
		boost::fusion::for_each(P, p);
	    }
	}

    };

    template<int N, class R>
    void evaluate(const Basis &basis, const R &r, vector<N> f) {
	f.clear_range(0, basis.size());
	foreach (const Basis::Block &block, basis.blocks()) {

#define DFT_WAVEFUNCTION_SHELL(R, DATA, TYPE)				\
	    if (block.shell().type == basis::shell::TYPE) {		\
		shell::evaluate<basis::shell::TYPE>(block, r, f);	\
		continue;						\
	    } else

	    BOOST_PP_SEQ_FOR_EACH(DFT_WAVEFUNCTION_SHELL,
				  (), BASIS_SHELL_TYPES);

#undef DFT_WAVEFUNCTION_SHELL


	}
    }

} // namespace wavefunction
} // namespace dft

#endif // DFT_WAVEFUNCTION_HPP
