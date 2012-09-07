#ifndef PARALLEL_HPP
#define PARALLEL_HPP

#include "config.h"

#ifdef HAVE_GA
#ifdef HAVE_MPI
#include <mpi.h>
#include "ga-mpi.h"
#endif
#endif


#include <iostream>
#include <memory>
#include <boost/noncopyable.hpp>
#include <boost/type_traits/is_scalar.hpp>
#include <boost/utility/enable_if.hpp>


namespace parallel {

#ifdef HAVE_MPI
#ifdef HAVE_GA
    inline MPI_Comm mpi_comm() {
	return GA_MPI_Comm();
    }
#endif
#endif

    struct Stream {
	explicit Stream(std::ostream *os = &std::cout)
	    : os_(os) {}
	operator bool() const { return (os_ != NULL); }
	Stream& operator<<(std::ostream& (*f)(std::ostream&)) {
	    if (os_) f(*os_);
	    return *this;
	}
	template<class O>
	Stream& operator<<(const O &o) {
	    if (os_) *os_ << o;
	    return *this;
	}
    private:
	std::ostream *os_;
    };

    struct Comm : boost::noncopyable {
	virtual ~Comm() {}
	virtual size_t rank() const = 0;
	virtual size_t size() const = 0;
	virtual void barrier() const = 0;
	virtual void reduce(std::string op, int *data, size_t size) const = 0;
	virtual void reduce(std::string op, double *data, size_t size) const = 0;
	virtual void broadcast(double *data, size_t size, int root) const = 0;
	Stream cout() const {
	    return Stream((rank() == 0) ? &std::cout : NULL);
	}
    };

    struct Task : boost::noncopyable {
	virtual ~Task() {}
	virtual size_t next() = 0;
	virtual void reset() = 0;
    };

    // inline
    // std::ostream& operator<<(std::ostream &ostream, const Environment &environment) {
    // 	ostream << "parallel::Environment(";
    // 	ostream << environment.rank() << ",";
    // 	ostream << environment.size();
    // 	ostream << ")";
    // 	return ostream;
    // }

}

struct Parallel : parallel::Comm {
    typedef parallel::Comm Comm;

    struct Task : boost::noncopyable {
	void reset() { impl_->reset(); }
	size_t operator++(int) { return impl_->next(); }
    private:
	friend class Parallel;
	Task(parallel::Task *task) : impl_(task) {}
	std::auto_ptr<parallel::Task> impl_;
    };

    Parallel();
    size_t rank() const;
    size_t size() const;
    void barrier() const;

    void broadcast(int *data, size_t size, int root) const;
    void broadcast(size_t *data, size_t size, int root) const;
    void broadcast(double *data, size_t size, int root) const;

    template<class V>
    typename boost::enable_if< boost::is_scalar<V> >::type
    broadcast(V &v, int root) const {
	broadcast(&v, 1, root);
    }

    template<class V>
    typename boost::disable_if< boost::is_scalar<V> >::type
    broadcast(V &v, int root) const {
	//std::cout << &v[0] << " " << v.size() << std::endl;
	broadcast(&v[0], v.size(), root);
    }

    void reduce(std::string op, int *data, size_t size) const;
    void reduce(std::string op, double *data, size_t size) const;

    template<class V>
    typename boost::enable_if< boost::is_scalar<V> >::type
    reduce(std::string op, V &v) const {
	reduce(op, &v, 1);
    }

    template<class V>
    typename boost::disable_if< boost::is_scalar<V> >::type
    reduce(std::string op, V &v) const {
	reduce(op, &v[0], v.size());
    }

    Task& task() const { return task_; }
    Comm& node() { return *node_; }
    Comm& cluster() { return *cluster_; }

private:
    std::auto_ptr<Comm> node_;
    std::auto_ptr<Comm> cluster_;
    mutable Task task_;
};

#endif // PARALLEL_HPP
