#ifndef OMP_HPP
#define OMP_HPP

#ifdef _OPENMP
#include <omp.h>
#endif

//#include <cstdlib>
#include <boost/noncopyable.hpp>
#include <boost/type_traits/remove_reference.hpp>

#ifndef _OPENMP
extern "C" {
    static int omp_get_max_threads() { return 1; }
    static int omp_get_num_threads() { return 1; }
    static int omp_get_thread_num() { return 0; }
    static void omp_set_nested(int) {}
    static void omp_set_num_threads(int) {}
}
#endif

namespace omp {

    struct thread {
	operator int() const {
	    return omp_get_thread_num();
	}
    };

    inline bool master() {
	return (thread() == 0);
    }

    inline size_t num_threads() {
	return omp_get_num_threads();
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

    template<typename T, typename V = T>
    struct task : boost::noncopyable {
	typedef V value_type;
	typedef typename boost::remove_reference<T>::type R;
	explicit task(const R &r = R()) : data_(r) {}
	explicit task(R &r) : data_(r) {}
	void reset(const value_type &value = value_type()) {
#pragma omp barrier
#pragma omp master
	    data_ = value;
#pragma omp barrier
	}
	//operator T&() const { return value_; }
	value_type operator++(int) {
	    value_type next;
#pragma omp critical(task)
	    {
		next = data_++;
	    }
	    return next;
	}
    private:
	T data_;
    };


    struct lock : boost::noncopyable {
#ifndef _OPENMP
	void set() {}
	void unset() {}
#else
    	lock() {
	    omp_init_lock(&lock_);
	}
	~lock() {
	    omp_destroy_lock(&lock_);
	}
	void set() {
	    omp_set_lock(&lock_);
	}
	void unset() {
	    omp_unset_lock(&lock_);
	}
    private:
	omp_lock_t lock_;
#endif
    };

    struct scoped_lock {
	scoped_lock(omp::lock &lock)
	    : lock_(lock) { lock_.set(); }
	~scoped_lock() { lock_.unset(); }
    private:
	omp::lock &lock_;
    };

}

#endif // OMP_HPP
