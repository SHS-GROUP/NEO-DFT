#ifndef OMP_HPP
#define OMP_HPP

#include <cstdlib>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace omp {

    struct thread {
	operator int() const {
#ifdef _OPENMP
	    return omp_get_thread_num();
#else
	    return 0;
#endif
	}
    };

    inline size_t num_threads() {
#ifdef _OPENMP
	return omp_get_num_threads();
#else
	return 1;
#endif
    }

    inline void set_num_threads(int n) {
#ifdef _OPENMP
	omp_set_num_threads(n);
#endif
    }

    inline void set_dynamic(bool f) {
#ifdef _OPENMP
	omp_set_dynamic(f);
#endif
    }

    inline void set_nested(bool f) {
#ifdef _OPENMP
	omp_set_nested(f);
#endif
    }

    inline void set_serial() {
	set_num_threads(1);
	set_dynamic(false);
	set_nested(false);
    }


    struct task {
	explicit task(int initial = 0) : value_(initial) {}
	int operator++(int) {
	    int next;
#pragma omp critical(task)
	    {
		next = value_++;
	    }
	    return next;
	}
    private:
	int value_;
    };


}

#endif // OMP_HPP
