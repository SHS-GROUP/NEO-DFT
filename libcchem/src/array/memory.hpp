#ifndef ARRAY_MEMORY_HPP
#define ARRAY_MEMORY_HPP

#include "array/array.hpp"

#include <string>
#include <vector>
#include <numeric>
#include <stdexcept>

#include <boost/typeof/typeof.hpp>

struct Array<>::Memory : Array<>::Factory {
    typedef Array<>::Type Type;
    struct Implementation;

    Array<>::Implementation*
    allocate(Type t, size_t N, const size_t *dims, const char *name) const;

    std::string str() const { return "Memory"; }

private:
    friend class Array<>::Memory::Implementation;
};


struct Array<>::Memory::Implementation : Array<>::Implementation {
private:
    Type type_;
    std::vector<size_t> dims_;
    std::string name_;
    Array<>::Memory factory_;
    std::vector<char> data_;
    std::vector<size_t> strides_;

public:

    Implementation(Type type, size_t N, const size_t *dims, const char *name,
		   const Memory &f)
	: type_(type),
	  dims_(dims, dims+N),
	  name_(name),
	  factory_(f)
    {
	size_t start[N];
	std::fill(start, start+N, 0);
	if (!check_range(N, start, &this->dims_[0]))
	    throw std::range_error("zero-sized dimension");

	size_t size = dims[0];
	strides_.push_back(1);
	for (size_t i = 1; i < N; ++i) {
	    size *= dims[i];
	    strides_.push_back(strides_[i-1]*dims[i-1]);
	}
	data_.resize(size*type.size);
    }

    const Array<>::Factory& factory() const { return factory_; }

    void put(const void *buffer, const size_t *start, const size_t *stop) {
	copy(start, stop, *this, Put(static_cast<const char*>(buffer)));
    }
    void get(void *buffer, const size_t *start, const size_t *stop) const {
	copy(start, stop, *this, Get(static_cast<char*>(buffer)));
    }
    void wait() {}
    void open() {}
    void close() {}

private:

    static
    bool check_range(size_t N, const size_t *start, const size_t *finish) {
	for (size_t i = 0; i < N; ++i) {
	    if (finish[i] <= start[i]) return false;
	}
	return true;
    }

    struct Put {
	const char* data_;
	explicit Put(const char* data) : data_(data) {}
	void operator()(char *begin, char *end) {
	    size_t size = (end-begin);
	    std::copy(data_, data_+size, begin);
	    data_ += size;
	}
    };

    struct Get {
	char* data_;
	explicit Get(char* data) : data_(data) {}
	void operator()(const char *begin, const char *end) {
	    size_t size = (end-begin);
	    std::copy(begin, end, data_);
	    data_ += size;
	}
    };

    template<class A, class F>
    static void copy(const size_t *start, const size_t *finish, A &a, F copy) {
	int N = a.dims_.size();

	if (!check_range(N, start, finish))
	    throw std::range_error("invalid index range");
	if (!check_range(N, start, finish))
	    throw std::range_error("index out of bounds");

	int M = -1;
	size_t dims[N], strides[N];
	size_t s[N], f[N], index[N];

	// collapse contiguos indices
	for (int i = 0; i < N; ++i) {
	    if ((i == 0) || ((finish[i] - start[i]) != a.dims_[i])) {
		++M;
		dims[M] = a.dims_[i];
		s[M] = start[i];
		f[M] = finish[i];
	    }
	    else {
		size_t n = (finish[i] - start[i]);
		dims[M] *= n;
		f[M] *= n;
	    }
	}
	M += 1;
	index[0] = s[0];
	strides[0] = 1;
	for (int i = 1; i < M; ++i) {
	    index[i] = s[i];
	    strides[i] = strides[i-1]*dims[i-1];
	}

	// copy contiguos blocks
	do {
	    size_t S = a.type_.size;
	    size_t n = (*f - *s);
	    std::ptrdiff_t diff = std::inner_product(index, index+M, strides, 0);
	    assert(S*(diff+n) <= a.data_.size());
	    BOOST_AUTO(data, (&a.data_[0])+S*diff);
	    copy(data, data+n*S);
	    *index += n;
	} while (next(M, s, f, index));
    }
    
    template <typename T, typename U>
    static bool next(size_t N, T start, T finish, U index) {
	for (size_t i = 0; i < N-1; ++i) {
	    ++(*index);
	    if (*index < *finish) return true;
	    *index = *start;
	    ++start;
	    ++finish;
	    ++index;
	}
	++(*index);
	if (*index < *finish) return true;
	return false;
    }

};

inline Array<>::Implementation*
Array<>::Memory::allocate(Type t, size_t N, const size_t *dims, const char *name) const {
    return new Implementation(t, N, dims, name, *this);
}


#endif // ARRAY_MEMORY_HPP
