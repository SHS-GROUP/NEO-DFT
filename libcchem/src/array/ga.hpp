#ifndef ARRAY_GA_HPP
#define ARRAY_GA_HPP

#include "array/array.hpp"

#include <cstring>
#include <string>
#include <vector>
#include <numeric>
#include <boost/typeof/typeof.hpp>

#include "ga.h"
#include "macommon.h"

struct Array<>::GA : Array<>::Factory {
    typedef Array<>::Type Type;
    struct Implementation;
    Array<>::Implementation*
    allocate(Type t, size_t N, const size_t *dims, const char *name) const;
    std::string str() const { return "GA"; }
    void set_chunk(size_t N, const size_t * dims) {
	chunk_ = std::vector<size_t>(dims, dims + N);
    }
    const size_t* get_chunk() const {
	return (chunk_.empty() ? NULL : &chunk_[0]);
    }
private:
    std::vector<size_t> chunk_;
};


struct Array<>::GA ::Implementation : Array<>::Implementation {
private:
    Type type_;
    std::vector<size_t> dims_;
    std::string name_;
    Array<>::GA  factory_;
    int handle_;

public:

    Implementation(Type type, size_t N, const size_t *dims, const char *name,
		   const GA &f)
	: type_(type), dims_(dims, dims+N), factory_(f)
    {

	// std::cout <<  "Array<>::GA ::Implementation("
	// 	  << N << ","
	// 	  << "(" << dims[0] << "),"
	// 	  << name
	// 	  << ")" << std::endl;

	std::vector<int> ga_dims(dims, dims+N);
	// GA C api assumes C-order
	std::reverse(ga_dims.begin(), ga_dims.end());

	size_t size = 1;
	for (size_t i = 0; i < N; ++i) {
	    size *= this->dims_[i];
	}

	// std::cout << "GA: "
	// 	  << limit << " "
	// 	  << GA_Memory_limited() << " "
	// 	  << GA_Memory_avail() << " "
	// 	  << size*type.size << std::endl;

	this->handle_ = 0;
	{ //if ((limit >= size*type.size) || !GA_Memory_limited()) {
	    char *ga_name = strdup(name);
	    std::vector<int> ga_chunk(ga_dims);
	    ga_chunk[0] = 1;
	    // std::cout << ga_chunk[0] << " " << ga_chunk[1] << " "
	    	      // << ga_chunk[2] << " " << ga_chunk[3] << std::endl;
	    this->handle_ =
		NGA_Create(ga_type(type), N, &ga_dims[0], ga_name, &ga_chunk[0]);
	    free(ga_name);
	}
	if (!this->handle_) throw Error(__FILE__ ": NGA_Create failed");

	GA_Zero(this->handle_);
    }

    ~Implementation() {
	GA_Destroy(this->handle_);
    }

    const Array<>::Factory& factory() const { return factory_; }

    void put(const void *buffer, const size_t *start, const size_t *stop) {
	apply(*this, &NGA_Put, start, stop, const_cast<void*>(buffer));
    }
    void get(void *buffer, const size_t *start, const size_t *stop) const {
	apply(*this, &NGA_Get, start, stop, buffer);
    }

    void wait() { GA_Sync(); }

    bool parallel() const { return true; }

private:
    static int ga_type(Type t) {
	if (t.id == Array<>::Type::DOUBLE) return C_DBL;
	if (t.id == Array<>::Type::INT) return C_INT;
	throw Error(__FILE__ ": unsupported GA type");
    }

    template<class A, class F>
    static size_t apply(A &ga, const F &f,
			const size_t *start, const size_t *finish,
			void *buffer) {

	size_t N = ga.dims_.size();

	// GA/armci uses int internally,
	// therefore operations over 2 gigabytes may overflow.
	// The logic below splits such operations.

	// split contiguous operations over 1GB
	size_t max = (1LU << 30)/ga.type_.size;
	size_t size = 1;
	for (size_t i = 0; i < N; ++i) {
	    size *= (finish[i] - start[i]);
	}

	// determine leading dimension
	std::vector<size_t> ld;
	ld.push_back(1);
	for (size_t i = 0; i < N; ++i) {
	    ld.push_back(ld[i]*ga.dims_[i]);
	}

	// find last non-one index
	size_t index = N;
	while (index) {
	    --index;
	    if ((finish[index] - start[index]) > 1) break;
	}

	// buffer is small enough
	if (size <= max && ld[index] <= max) 
	    return apply(f, ga.handle_, N, start, finish, buffer);

	bool contiguous = true;
	for (size_t i = 0; i <= index; ++i) {
	    if (finish[i] - start[i] != ga.dims_[i]) contiguous = false;
	}

	size_t chunk = 1;
	// entire buffer is contiguous
	if (contiguous) {
	    size_t n = size/(finish[index] - start[index]);
	    chunk = std::max<size_t>(max/n, 1);
	}

	size_t count = 0;
	std::vector<size_t> s_i(start, start+N); 
	std::vector<size_t> f_i(finish, finish+N); 
	for (size_t i = start[index]; i < finish[index]; i += chunk) {
	    s_i[index] = i;
	    f_i[index] = std::min(i+chunk, finish[index]);
	    char* ptr = ((char*)buffer) + count*ga.type_.size;
	    count += apply(ga, f, &s_i[0], &f_i[0], ptr);
	}
	CCHEM_ASSERT(count == size);
	return count;
    }

    template<class F>
    static size_t apply(const F &f, int ga, size_t N,
			const size_t *start, const size_t *finish,
			void *buffer) {
	size_t size = 1;
	std::vector<int> lo, hi, ld;
	// GA C api assumes C-order
	start += N;
	finish += N;
	for (int i = N; i > 0; --i) {
	    --start;
	    --finish;
	    lo.push_back((*start));
	    hi.push_back((*finish)-1);
	    ld.push_back(*finish - *start);
	    size *= ld.back();
	}
#pragma omp critical(comm)
	f(ga, &lo[0], &hi[0], buffer, &ld[1]);
	return size;
	
    }

};

inline Array<>::Implementation*
Array<>::GA ::allocate(Type t, size_t N, const size_t *dims, const char *name) const {
    return new Implementation(t, N, dims, name, *this);
}


#endif // ARRAY_GA_HPP
