#ifndef ARRAY_ARRAY_HPP
#define ARRAY_ARRAY_HPP

#include <string>
#include <vector>

#include <boost/mpl/vector.hpp>
#include <boost/mpl/find.hpp>
#include <boost/mpl/at.hpp>

#include <boost/array.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/shared_ptr.hpp>



#include "boost/utility/timer.hpp"


template<typename T = void>
struct Array;

template<>
class Array<void> {

public:

    template<typename T = void>
    struct Type {
	typedef boost::mpl::vector<
	    void,
	    double
	    > vector;

	static const int value = boost::mpl::distance<
	    typename boost::mpl::begin<vector>::type,
	    typename boost::mpl::find<vector, T>::type
	    >::value;
    };

    enum TYPE {
	VOID = Type<void>::value,
	DOUBLE = Type<double>::value
    };


private:
    struct Implementation;

public:
    struct GA;
    struct HDF5;

    struct Factory {
	virtual
	Implementation* allocate(TYPE type,
				 size_t N, const size_t *dims,
				 const char *name) const = 0;
    };    

private:
    struct Implementation : boost::noncopyable {
	virtual
	void put(const void *buffer, const size_t *start, const size_t *stop) = 0;
	virtual
	void get(void *buffer, const size_t *start, const size_t *stop) const = 0;
	virtual void synchronize() = 0;
	virtual const Factory& factory() const = 0;
    };

public:
    struct Counters {
	typedef boost::utility::timer timer;
	typedef timer::value_type value_type;
	void clear() { *this = Counters(); }
	const value_type& read() const { return read_; }
	const value_type& write() const { return write_; }
    private:
	template<typename> friend class Array;
	value_type read_, write_;
    };

};


template<typename T>
struct Array : boost::noncopyable {

    typedef size_t size_type;
    typedef Array<>::Factory Factory;

    const size_type N;

    template<typename S>
    Array(size_type N, const S &dims, const char *name, const Factory &f)
	: N(N), dims(&dims[0], &dims[0]+N)
    {
	Array<>::TYPE type = Array<>::TYPE(Array<>::Type<T>::value);
	size_type r[N];
	copy(N, dims, r);
	impl_.reset(f.allocate(type, N, r, name));
    }
    virtual ~Array() {}

    const size_type* shape() const { return &dims[0]; }

    const Factory& factory() const { return impl_->factory(); }

    template<typename S>
    void put(const T *buffer, const S &start, const S &stop) {
	size_type r1[N], r2[N];
	copy(N, start, r1);
	copy(N, stop, r2);
	impl_->put(buffer, r1, r2);
    }

    template<typename S>
    void get(T *buffer, const S &start, const S &stop) const {
	size_type r1[N], r2[N];
	copy(N, start, r1);
	copy(N, stop, r2);
	impl_->get(buffer, r1, r2);
    }

    struct invalid_dimension {};

    template<typename S>
    T at(const S &index) const {
	size_type start[N], stop[N];
	for (size_t i = 0; i < N; ++i) {
	    start[i] = index[i];
	    stop[i] = index[i]+1;
	}
	T value;
	impl_->get(&value, start, stop);
	return value;
    }

    T at(int i, int j, int k, int l) const {
	if (N != 4) throw invalid_dimension();
	size_type index[] = { i,j,k,l };
	return at(index);
    }

    void synchronize() {
	impl_->synchronize();
    }

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

protected:
    std::vector<size_type> dims;
    std::auto_ptr<Array<>::Implementation> impl_;

private:
    template<typename S>
    static void copy(size_type N, const S &it, size_type *ot) {
	std::copy(&it[0], &it[0]+N, ot);
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



#endif // ARRAY_ARRAY_HPP
