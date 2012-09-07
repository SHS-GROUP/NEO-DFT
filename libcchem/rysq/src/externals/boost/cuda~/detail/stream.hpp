#ifndef BOOST_CUDA_DETAIL_STREAM_HPP
#define BOOST_CUDA_DETAIL_STREAM_HPP

#include "boost/cuda/stream.hpp"
#include "boost/cuda/exception.hpp"
#include <boost/type_traits/remove_pointer.hpp>
#include <boost/shared_ptr.hpp>

#include <cuda_runtime.h>


namespace boost {
namespace cuda {

    template<>
    cudaStream_t stream::data<cudaStream_t>() const {
	return cudaStream_t(data());
    }

    stream::type stream::create_() {
	cudaStream_t s;
	cudaStreamCreate(&s);
	check_status();//cuda_throw_error( );
	//std::cout << "create " << s << std::endl;
	return type(s);
    }

    void stream::destroy_(type data) {
	if (data) {
	    //std::cout << "destroy " << data << std::endl;
	    cudaStreamDestroy(cudaStream_t(data));
	    check_status();//cuda_throw_error( );
	}
    }

    struct stream::impl {
	//typedef boost::remove_pointer<cudaStream_t>::type type;
	boost::shared_ptr<type> data_;
	impl(type data)
	    : data_(new type(static_cast<type>(data)), destructor()) {}	
	type get() const { return *data_.get(); }
	struct destructor {
	    void operator()(type *data) const {
		stream::destroy_(*data);
		delete data;
	    }
	};
    };   

    stream::stream() {
	data_ = new impl(create_());
	synchronous_ = 0;
    }

    stream::stream(type data) {
	data_ = (!data) ? 0 : new impl(data);
	synchronous_ = 0;
    }

    stream::~stream() {
	if (data_) {
	    if (data_) delete data_;
	}
    }

    void stream::synchronize() {
	cudaStreamSynchronize(this->data<cudaStream_t>());
	check_status();
    }

    void stream::copy(const stream &s) {
	impl *old = data_;
	if (s.data_) {
	    data_ = new impl(*s.data_);
	}
	if (old) { delete old; }
	synchronous_ = s.synchronous_;    
    }    

    stream::type stream::data() const {
	if (!data_) return 0;
	return data_->get();
    }


}
}

#endif // BOOST_CUDA_DETAIL_STREAM_HPP
