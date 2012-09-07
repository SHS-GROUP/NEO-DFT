#ifndef CCHEM_HF_SCREEN_HPP
#define CCHEM_HF_SCREEN_HPP

#include "hf/hf.hpp"

#include "boost/utility/profiler.hpp"
#include <string>

namespace hf {

    template<class M1, class M2>
    void Dmax(const Basis &basis, const M1 &D, M2 &Dmax) {
	for (int b = 0; b < basis.num_shells(); ++b) {
	    int nj = basis.shell(b).size();
	    int j0 = basis.shell(b).start();//.firstShellFunction(b);

	    for (int a = 0; a <= b; ++a) {
		int ni = basis.shell(a).size ();
		int i0 = basis.shell(a).start();//.firstShellFunction(a);

		double q = 0.0;
		for (int j = j0; j < j0+nj; ++j) {
		    for (int i = i0; i < i0+ni; ++i) {
			q = std::max(q, fabs(D(i,j)));
		    }
		}
		Dmax(a,b) = q;
		Dmax(b,a) = q;

	    }
	}
    }

    template<class M1, class M2 = M1>
    class Screen {
    public:
	Screen(const M1 &Kmax, const M2 &Dmax, double cutoff) :
	    Kmax_(Kmax), Dmax_(Dmax), cutoff_(cutoff) {}

	template<typename T>
	static T max(const T &a, const T &b, const T &c, const T &d) {
	    return std::max(std::max(a,b), std::max(c,d));
	}
	template<class M>
	static double max(const M &Dmax, int i, int j, int k, int l) {
	    return std::max((4.0*std::max(Dmax(i,j), Dmax(k,l))),
			    (max(Dmax(i,k), Dmax(i,l), Dmax(j,k), Dmax(j,l))));
	}

	struct Functor {
	private:
	    const M1 &Kmax_;
	    const M2 &Dmax_;
	    double cutoff_;
	    size_t max_;
	    int begin_[4], end_[4];
	    mutable int restart_[4];
	    bool ab_, ac_, cd_;
	    size_t size_;
	    std::string name_;
	public:
	    Functor(const M1 &Kmax, const M2 &Dmax, double cutoff, size_t max,
		    const Basis::Block &a, const Basis::Block &b,
		    const Basis::Block &c, const Basis::Block &d)
		: Kmax_(Kmax), Dmax_(Dmax), cutoff_(cutoff), max_(max)
	    {
		int begin[] = { a.start(), b.start(), c.start(), d.start() };
		int end[] = { a.stop(), b.stop(), c.stop(), d.stop() };
		for (int i = 0; i < 4; ++i) {
		    this->begin_[i] = begin[i];
		    this->end_[i] = end[i];
		    this->restart_[i] = begin[i];
		}
		ab_ = (a == b);
		ac_ = (a == c);
		cd_ = (c == d);

		size_ = (a.shell().size()*b.shell().size()*
			 c.shell().size()*d.shell().size());
		// BOOST_PROFILE_LINE;
		// std::ostringstream os;
		// os << "(" << a << b << "|" << c << d << ")";
		// name_ = (os.str() == "(ss|ss)") ? "(ss|ss)" : "(xx|xx)";
	    }
	    double cutoff() const { return cutoff_; }
	    template<class C, class F>
	    bool operator()(C &quartets, const F &f) const {
		return apply<false>(quartets, f);
	    }
	    template<class C>
	    bool operator()(C &quartets) const {
		return apply<false>(quartets, Identity());
	    }
	private:
	    friend class Screen;

	    void update_profiler(size_t size) const {
		//boost::utility::global_profiler()[name_] += size*size_;
	    }

	    struct Identity {
		template<class C>
		C operator()(const C &c) const { return c; }
	    };

	    template<bool test_only, class C, class F>
	    bool apply(C &quartets, const F &f) const {
		BOOST_PROFILE_LINE;
		// build vector of quartets
		for (int l = restart_[3]; l < end_[3]; ++l) {
		    for (int j = restart_[1]; j < end_[1]; ++j) {
			restart_[1] = begin_[1];

			int kfirst = (cd_) ? l : begin_[2];
			for (int k = kfirst; k < end_[2]; ++k) {
			    int ifirst = (ab_) ? j : begin_[0];
			    int ilast = (ac_) ? k + (j <= l) : end_[0];
			    for (int i = ifirst; i < ilast; ++i) {
				if (!test(i, j, k, l)) continue;
				if (test_only) return true;
				typename C::value_type quartet =
				    {{i, j, k, l}};
				quartets.push_back(f(quartet));
			    }
			}
			if (quartets.size() > max_) {
			    restart_[3] = l;
			    restart_[1] = j+1;
			    if (restart_[1] == end_[1]) {
				restart_[1] = begin_[1];
				restart_[3]++;
			    }
			    update_profiler(quartets.size());
			    return true;
			}
		    }
		}
		//update_profiler(quartets.size());
		restart_[3] = end_[3];
		return ((test_only) ? false : quartets.empty());
	    }

	    bool test(int i, int j, int k, int l) const {
		int f = 1;
		return (f*max(Dmax_,i, j, k, l)*Kmax_(i,j)*Kmax_(k,l) > cutoff_);
	    }
	};

	template<class Task>
	Functor operator()(const std::vector<Basis::Block> &blocks,
			   const Task &task,
			   size_t max  = std::numeric_limits<size_t>::max()) const {
	    return (*this)(blocks.at(task[0]), blocks.at(task[1]),
			   blocks.at(task[2]), blocks.at(task[3]), max);
	}

	Functor operator()(const Basis::Block &a, const Basis::Block &b,
			   const Basis::Block &c, const Basis::Block &d,
			   size_t max  = std::numeric_limits<size_t>::max()) const {
	    // BOOST_PROFILE_LINE;
	    return Functor(Kmax_, Dmax_, cutoff_, max, a,b,c,d);
	}

	template<typename  iterator>
	bool test(const iterator (&block)[4]) const {
	    return test(*block[0], *block[1], *block[2], *block[3]);
	}

	bool test(const Basis::Block &a, const Basis::Block &b,
			const Basis::Block &c, const Basis::Block &d)  const {
	    std::vector<boost::array<bool,4> > quartets;
	    return this->operator()(a,b,c,d).apply<true>(quartets);
	}

	double cutoff() const { return cutoff_; }
    private:
	const M1 &Kmax_;
	const M2 &Dmax_;
	double cutoff_;
    };

}

#endif /* CCHEM_HF_SCREEN_HPP */
