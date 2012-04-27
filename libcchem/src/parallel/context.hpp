#ifndef PARALLEL_CONTEXT_HPP
#define PARALLEL_CONTEXT_HPP

#include <iostream>

namespace parallel {

    struct Context {
	struct tuple {
	    explicit tuple(int rank = 0, int size = 1)
		: rank_(rank), size_(size) {}
	    int rank() const { return rank_; }
	    int size() const { return size_; }
	private:
	    int rank_, size_;
	};
	typedef tuple Node;
	typedef tuple SMP;
	Context();
	int rank() const { return rank_; }
	int size() const { return size_; }
	const SMP& smp() const { return smp_; }
    private:
	int rank_, size_;
	SMP smp_;
    };

    inline
    std::ostream& operator<<(std::ostream &ostream, const Context &context) {
	ostream << "parallel::Context(";
	ostream << context.rank() << ",";
	ostream << context.size();
	ostream << ")";
	return ostream;
    }

}

#endif // PARALLEL_CONTEXT_HPP
