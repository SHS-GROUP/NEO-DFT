#ifndef CC_TRIPLES_DEVICE_HPP
#define CC_TRIPLES_DEVICE_HPP

#include "cc/triples/triples.hpp"
#include "cc/triples/thread.hpp"

#include <cuda.h>
#include "blas.hpp"

namespace cc {
namespace triples {
namespace detail {

    struct Device
	: detail::Thread< array::adapter<double, 3, array::device_tag> >
    {
      typedef array::adapter<double, 3, array::device_tag> T3;
	typedef detail::Thread<T3> Thread;
      
    private:
	Matrix tmp_, hvvv_;
	cublas::matrix<double> t2_, Vjk_, Vkj_, vov_, vvv_;
	cublas::matrix<double> v_[3];
	std::vector<cudaStream_t> streams_;
	double *H_;

	static size_t num_streams(int device) {
	    cudaDeviceProp p;
	    cudaGetDeviceProperties(&p, device);
	    return 8; //p.deviceOverlap;
	}

	static void* allocate_impl(size_t size);
	static void free_impl(void* ptr);

    public:
	template<size_t N, typename U>
	static array::adapter<double, N, array::device_tag>
	make_array(const U &shape, double *data = 0) {
	    size_t size = 1;
	    for (size_t i = 0; i < N; ++i)
		size *= shape[i];
	    if (!data) 
	      data = Device::template allocate<double>(size);
	    return array::adapter<double, N, array::device_tag>(data, shape);
	}
	template<typename T>
	static T* allocate(size_t size) {
	    return (T*)(allocate_impl(size*sizeof(T)));
	}
	template<typename T>
	static void free(T *ptr) {
	    free_impl(ptr);
	}

    public:
	size_t done;

	Device(const Thread &thread, size_t ns = 0);

	~Device();

	void ijk1_impl(const Data &data,
		       const Matrix &tab, const Matrix &taib,
		       const Matrix &Vjk, const Matrix &Vkj,
		       T3 &t3,
		       int i, int j, int k,
		       std::pair<bool,int> terms,
		       Queue &queue);
    };
}
}
}

#endif // CC_TRIPLES_DEVICE_HPP
