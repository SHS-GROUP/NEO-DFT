#include <rysq-cuda.hpp>

#include "cuda.hpp"
#include <assert.h>
#include <map>

#include <boost/current_function.hpp>
#include <boost/typeof/typeof.hpp>
#include "foreach.hpp"

namespace rysq {
namespace cuda {
namespace runtime {

    struct Index {
	typedef boost::array<rysq::type,4> key;
	Index() : data_(0) {}
	~Index() { reset(); }
	void reset() {
	    if (data_) ::cuda::free(data_);
	    data_ = 0;
	}
	const ushort4* find(rysq::type A, rysq::type B,
			    rysq::type C, rysq::type D) {
	    if (!data_) generate();
	    key key = {{ A, B, C, D }};
	    BOOST_AUTO(it, map_.find(key));
	    if (it == map_.end()) return NULL;
	    return data_ + it->second;
	}
    private:
	ushort4 *data_; 
	std::map<key,size_t> map_;
	void generate() {
	    size_t index = 0;
	    std::vector<ushort4> host;
	    foreach (rysq::type d, rysq::types) {
		foreach (rysq::type c, rysq::types) {
		    foreach (rysq::type b, rysq::types) {
			foreach (rysq::type a, rysq::types) {
			    rysq::shell A(a);
			    rysq::shell B(b);
			    rysq::shell C(c);
			    rysq::shell D(d);
			    size_t size = A.size*B.size*C.size*D.size;
			    if (A.L+B.L+C.L+D.L > 9) continue;
			    append(A, B, C, D, host);
			    key key = {{ a, b, c, d }};
			    //std::cout << size << std::endl;
			    map_[key] = index;
			    index += size;
			}
		    }
		}
	    }
	    data_ = ::cuda::malloc<ushort4>(host.size());
	    ::cuda::copy(host.size(), &host[0], data_,
			 ::cuda::host_to_device);
	}
	void append(rysq::shell a, rysq::shell b,
		    rysq::shell c, rysq::shell d,
		    std::vector<ushort4> &data) {
	    int m1 = (a.L+1);
	    int m2 = (b.L+1)*m1;
	    int m3 = (c.L+1)*m2;
	    for (int l = d.begin(); l < d.end(); ++l) {
	    for (int k = c.begin(); k < c.end(); ++k) {
	    for (int j = b.begin(); j < b.end(); ++j) {
	    for (int i = a.begin(); i < a.end(); ++i) {
		ushort4 idx = {};
		idx.x = LX[i] + m1*LX[j] + m2*LX[k] + m3*LX[l];
		idx.y = LY[i] + m1*LY[j] + m2*LY[k] + m3*LY[l];
		idx.z = LZ[i] + m1*LZ[j] + m2*LZ[k] + m3*LZ[l];
		// sp mask
		idx.w = 0;
		idx.w |= ((a == rysq::SP) && (i == 0)) << 0;
		idx.w |= ((b == rysq::SP) && (j == 0)) << 1;
		idx.w |= ((c == rysq::SP) && (k == 0)) << 2;
		idx.w |= ((d == rysq::SP) && (l == 0)) << 3;
		data.push_back(idx);
	    }
	    }
	    }
	    }
	}
    };



    void* malloc(size_t size) {
	void *ptr = ::cuda::malloc<char>(size);
	// std::cout << BOOST_CURRENT_FUNCTION << ": "
	// 	  << ptr << " "
	// 	  << size << " bytes" << std::endl;
	return ptr;
    }

    void free(void *ptr) {
	//std::cout <<  BOOST_CURRENT_FUNCTION << ": " << ptr << std::endl;
	::cuda::free(ptr);
    }

    void memset(void *ptr, int value, size_t size) {
	::cuda::memset(ptr, value, size);
    }

    template<>
    void copy(size_t size, const void *from, void *to, copy_type t) {
	// std::cout <<  BOOST_CURRENT_FUNCTION << ": "
	// 	  << from << " to " << to << " "
	// 	  << size << " bytes " << std::endl;
	if (t == device_to_host)
	    ::cuda::copy(size, from, to, ::cuda::device_to_host);
	else if (t == device_to_device)
	    ::cuda::copy(size, from, to, ::cuda::device_to_device);
	else if (t == host_to_device)
	    ::cuda::copy(size, from, to, ::cuda::host_to_device);
	else {
	    throw std::runtime_error(BOOST_CURRENT_FUNCTION);
	}
    }


    template<typename T>
    void copy(size_t size, const T *from, T *to, copy_type t,
	      const Stream &s);

    template<>
    void copy(size_t size, const void *from, void *to, copy_type t,
	      const Stream &stream) {
	cudaStream_t s = reinterpret_cast<cudaStream_t>(stream.data());
	if (t == device_to_host)
	    ::cuda::copy(size, from, to, ::cuda::device_to_host, s);
	else if (t == device_to_device)
	    ::cuda::copy(size, from, to, ::cuda::device_to_device, s);
	else if (t == host_to_device)
	    ::cuda::copy(size, from, to, ::cuda::host_to_device, s);
	else {
	    throw std::runtime_error(BOOST_CURRENT_FUNCTION);
	}
    }


    void Stream::synchronize() {
	::cuda::stream(static_cast<cudaStream_t>(data_)).synchronize();
    }

    Context::Context() {
	index_.reset(new Index());
    }

    Context::~Context() {}

    const cudaDeviceProp& Context::properties() const {
	if (!properties_.get()) {
	    properties_.reset(new cudaDeviceProp());
	    int device;
	    CUDA_VERIFY(cudaGetDevice(&device));
	    CUDA_VERIFY(cudaGetDeviceProperties(properties_.get(), device));
	}
	return *properties_;
    }

    const void* Context::index(rysq::type A, rysq::type B,
			       rysq::type C, rysq::type D) const {
	return index_->find(A,B,C,D);
    }
    
}
}
}

