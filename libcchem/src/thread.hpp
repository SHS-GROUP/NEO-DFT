#ifndef CCHEM_THREAD_HPP
#define CCHEM_THREAD_HPP

#include "runtime.hpp"
#include "blas.hpp"
#include "config.h"
#ifdef HAVE_CUBLAS
#include "cublas.hpp"
#endif

#include "foreach.hpp"
#include "utility/timer.hpp"

#include "cc/tensor.hpp"
#include "cc/utility.hpp"

#include <vector>
#include <map>
#include <set>
#include <memory>
#include <boost/noncopyable.hpp>

namespace cchem {

#ifdef HAVE_CUBLAS
    struct Device {
	cublas::handle_t cublas_handle;
	typedef cuda::event Event;
	struct Stream : cuda::stream {
	    cc::Symbol< cc::tensor_reference<4> > S;
	    cc::Symbol< cc::tensor_reference<4> > H;
	    Event event;
	};
	std::vector<Stream> streams;
	std::vector<Event> events;
	cc::Symbol< cc::tensor_reference<4> > S;
	cc::Symbol< cc::tensor_reference<4> > H;
	explicit Device(size_t ns = 0) {
	    this->cublas_handle = cublas::create();
	    stream_ = 0;
	    streams.resize(ns);
	    events.resize(ns);
	    typedef std::pair<const int,double*> P;
	    foreach (Stream &s, streams) {
		s.create();
		s.event.create();
	    }
	    foreach (Event &e, events) {
		e.create();
	    }
	}
	~Device() {
	    this->free();
	    cublas::destroy(this->cublas_handle);
	    foreach (void* ptr, host_) {
		cuda::free_host(ptr);
	    }
	    foreach (Stream &s, streams) {
		s.destroy();
		s.event.destroy();
	    }
	    foreach (Event &e, events) {
		e.destroy();
	    }
	}
	/** free device memory */
	void free() {
	    foreach (void *ptr, this->device_) {
		cuda::free(ptr);
	    }
	    device_.clear();
	}
	double* malloc(size_t size) {
	    double *ptr = cuda::malloc<double>(size);
	    device_.insert(ptr);
	    return ptr;
	}
	double* malloc_host(size_t size) {
	    double *ptr = cuda::malloc_host<double>(size);
	    host_.insert(ptr);
	    return ptr;
	}
	Stream& stream() {
	    return streams.at(stream_%streams.size());
	}
	void next() {
	    for (size_t i = 0; i < streams.size(); ++i) {
		if (streams.at(i).query()) {
		    stream_ = i;
		    return;
		}
	    }
	    ++stream_;
	}
    private:
	std::set<void*> host_;
	std::set<void*> device_;
	size_t stream_;
    };
#endif // HAVE_CUBLAS


    struct Thread : boost::noncopyable {
	template<class Counter>
	struct Task;

	cc::Buffer<double> buffer;
	cc::Symbol< cc::tensor_reference<4> > S;
	std::map<int,double*> data;

	static double* realloc(double* ptr, size_t size) {
	    return Runtime::rt().memory().realloc<double>(ptr, size);
	}
	static double* malloc(size_t size) {
	    return Runtime::rt().memory().malloc<double>(size);
	}
	static void free(double *ptr) {
	    Runtime::rt().memory().free(ptr);
	}
	void free() {
	    BOOST_AUTO(it, data.begin());
	    while (it != data.end()) {
		Runtime::rt().memory().free(it->second);
		++it;
	    }
	    data.clear();
	}

	Thread() : device(NULL) {}
	~Thread() {
	    delete device;
	    this->free();
	}

#ifndef HAVE_CUBLAS
	struct Device {};
#endif // HAVE_CUBLAS
	Device* device;

    };


    template<class Counter>
    struct Thread::Task {
    private:
	Counter counter_;
	int value_;
    public:
	int index;
	template<class C>
	explicit Task(C &c) : counter_(c) {}
	operator int() const { return value_; }
	int operator++(int) {
	    return next();
	}
	int operator++() {
	    return next();
	}
	void reset() {
#pragma omp barrier
#pragma omp master
	    counter_.reset();
#pragma omp barrier
	    value_ = 0;
	    index = 0;
	}
    private:
	int next() {
	    value_ = counter_++;
	    return value_;
	}
    };

} // namespace

#endif // CCHEM_THREAD_HPP

