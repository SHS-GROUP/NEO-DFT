#ifndef ARRAY_ARRAY_HPP
#define ARRAY_ARRAY_HPP

#include <string>
#include <vector>

#include <boost/mpl/vector.hpp>
#include <boost/mpl/find.hpp>
#include <boost/mpl/at.hpp>

#include <boost/type.hpp>
#include <boost/array.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/shared_ptr.hpp>

#include "parallel.hpp"
#include "utility/timer.hpp"


template<typename T = void>
struct Array;

template<>
class Array<void> {

public:

    struct Type {

	typedef boost::mpl::vector<
	    void,
	    int,
	    double
	    > vector;

	template<typename T>
	struct find {
	    static const int value = boost::mpl::distance<
		typename boost::mpl::begin<Type::vector>::type,
		typename boost::mpl::find<Type::vector, T>::type
		>::value;
	};

	enum TYPE {
	    VOID = Type::find<void>::value,
	    INT = Type::find<int>::value,
	    DOUBLE = Type::find<double>::value
	};

	const int id;
	const size_t size;
	template<typename T>
	Type(boost::type<T>) : id(find<T>::value), size(sizeof(T)) {}
    };

    struct Error : std::runtime_error {
	Error(const char *what) : std::runtime_error(what) {}
    };

    struct Counters {
	typedef utility::timer timer;
	timer::value_type read() const { return read_; }
	timer::value_type write() const { return write_; }
	void reset() {
	    read_ = timer::value_type();
	    write_ = timer::value_type();
	}
    private:
	template<typename> friend class Array;
	utility::timer::value_type read_, write_;
    };

private:
    struct Implementation;

public:
    struct Memory;
    struct GA;
    struct HDF5;

    struct Factory {
	virtual ~Factory() {}
	virtual
	Implementation* allocate(Type type,
				 size_t N, const size_t *dims,
				 const char *name) const = 0;
	virtual std::string str() const = 0;
    };

private:
    struct Implementation : boost::noncopyable {
	virtual ~Implementation() {}
	virtual
	void put(const void *buffer, const size_t *start, const size_t *stop) = 0;
	virtual
	void get(void *buffer, const size_t *start, const size_t *stop) const = 0;
	virtual
	void get_all(void *buffer, const size_t *start, const size_t *stop) const {
	    get(buffer, start, stop);
	}
	virtual void close() const {}
	virtual void open() const {}
	virtual void wait() = 0;
	virtual bool parallel() const { return false; }
    };

};


template<typename T>
struct Array : boost::noncopyable {

public:

    typedef size_t size_type;

    typedef Array<>::Type Type;
    typedef Array<>::Counters Counters;
    typedef Array<>::Factory Factory;
    typedef Array<>::Error Error;

    typedef Array<>::Memory Memory;
    typedef Array<>::HDF5 HDF5;
    typedef Array<>::GA GA;

    const size_type N;

    template<typename S>
    Array(size_type N, const S &dims, const char *name, const Factory &f)
	: N(N), dims(&dims[0], &dims[0]+N)
    {
	size_type r[N];
	copy(N, dims, r);
	impl_.reset(f.allocate(Type(boost::type<T>()), N, r, name));
	str_ = ("Array::" + f.str() + ": " + name);
    }
    virtual ~Array() {}

    const size_type* shape() const { return &dims[0]; }

    template<typename S>
    void put(const T *buffer, const S &start, const S &stop) {
	size_type r1[N], r2[N];
	copy(N, start, r1);
	copy(N, stop, r2);
	Counters::timer t;
	impl_->put(buffer, r1, r2);
	counters_.write_ += t;
    }

    template<typename S>
    void get(T *buffer, const S &start, const S &stop) const {
	get(buffer, start, stop, false);
    }

    template<typename S>
    void get_all(T *buffer, const S &start, const S &stop) const {
	get(buffer, start, stop, true);
    }

    // void put(const T *buffer, const size_t *start, const size_t *stop) {
    // 	put<size_t>(buffer, start, stop);
    // }

    // void get(T *buffer, const size_t *start, const size_t *stop) const {
    // 	get<size_t>(buffer, start, stop);
    // }

    struct invalid_dimension {};

    // template<typename R>
    // void get(T *buffer, const R &r1, const R &r2, const R &r3, const R &r4) {
    // 	if (N != 4) throw invalid_dimension();
    // 	size_type start[] = { r1.start(), r2.start(), r3.start(), r4.start() };
    // 	size_type stop[] = { r1.size(), r2.size(), r3.size(), r4.size() };
    // 	for (size_t i = 0; i < 4; ++i) {
    // 	    stop[i] += start[i];
    // 	}
    // 	get(buffer, start, stop);
    // }

    template<typename S>
    T at(const S &index) const {
	size_t start[N], stop[N];
	for (size_t i = 0; i < N; ++i) {
	    start[i] = index[i];
	    stop[i] = index[i]+1;
	}
	T value;
	get<size_t*>(&value, start, stop);
	return value;
    }

    T at(int i, int j, int k, int l) const {
	if (N != 4) throw invalid_dimension();
	size_type index[] = { i,j,k,l };
	return at(index);
    }

    void flush() {
	impl_->wait();
    }

    bool parallel() const {
	return impl_->parallel();
    }

    Counters& counters() const { return counters_; }

    const std::string& str() const {
	return str_;
    }

protected:
    std::vector<size_type> dims;
    std::auto_ptr<Array<>::Implementation> impl_;
    mutable Array<>::Counters counters_;
    std::string str_;

public:

    static boost::array<int,2> index(int i, int j) {
	boost::array<int,2> index = {{ i,j }};
	return index;
    }

    static boost::array<int,3> index(int i, int j, int k) {
	boost::array<int,3> index = {{ i,j,k }};
	return index;
    }

    static boost::array<int,4> index(int i, int j, int k, int l) {
	boost::array<int,4> index = {{ i,j,k,l }};
	return index;
    }

private:

    template<typename S>
    static void copy(size_type N, const S &it, size_type *ot) {
	std::copy(&it[0], &it[0]+N, ot);
    }

    template<typename S>
    void get(T *buffer, const S &start, const S &stop, bool all) const {
	size_type r1[N], r2[N];
	copy(N, start, r1);
	copy(N, stop, r2);
	Counters::timer t;
	if (all) impl_->get_all(buffer, r1, r2);
	else impl_->get(buffer, r1, r2);
	counters_.read_ += t;
    }

    void verify(const T *buffer, size_type *r1, size_type *r2) {
	std::cout << buffer << std::endl;
	for (size_type j = r1[3]; j < r2[3]; ++j) {
	    for (size_type i = r1[2]; i < r2[2]; ++i) {
		for (size_type b = r1[1]; b < r2[1]; ++b) {
		    for (size_type a = r1[0]; a< r2[0]; ++a) {
			size_type index[] = { a,b,i,j };
			T value = this->at(index);
			printf("verify %i %i %i %i: %13.10e == %13.10e\n",
			       int(a+1),int(b+1),int(i+1),int(j+1),
			       *buffer, value);
			assert(*buffer++ == value);
		    }
		}
	    }
	}
    }

};

template<typename T>
std::ostream& operator<<(std::ostream& os, const Array<T> &a) {
    os << a.str();
    size_t size = 1;
    os << " { ";
    for (size_t i = 0; i < a.N; ++i) {
	if (i > 0) os << ", ";
	os << a.shape()[i];
	size *= a.shape()[i];
    }
    os << " }, " << (sizeof(T)*size)/1e6 << " MB";
    os << " parallel=" << a.parallel();
    return os;
}

#endif // ARRAY_ARRAY_HPP
