#ifndef RYSQ_CUDA_DETAIL_HPP
#define RYSQ_CUDA_DETAIL_HPP


#include "rysq-eri.hpp"
#include "rysq-cuda.hpp"
#include "transpose.hpp"


#include <memory>
#include <boost/config.hpp>

#include <cuda.h>
#include <cuda_runtime.h>


namespace rysq {
namespace cuda {
namespace detail {


    struct Quartets {
	typedef const Int4* const_pointer;
	const Int4 *data_;
	size_t size_;
	Quartets(const cuda::Quartets &quartets)
	    : data_(quartets.begin()), size_(quartets.size()) {}
	const_pointer data()const { return data_; }
	BOOST_GPU_ENABLED
	size_t size() const { return size_; }
	BOOST_GPU_ENABLED 
	const Int4& operator[](int index) const {
	    return data_[index];
	}
    };

    struct Centers {
	typedef const rysq::Center* pointer;
	Centers() : data_() {}
	explicit Centers(const cuda::Centers &centers)
	    : data_(centers.begin()) {}
	BOOST_GPU_ENABLED
	operator pointer() const { return data_; }
    private:
	pointer data_;
    };


    struct Eri {
    	typedef rysq::Eri::Parameters Parameters;
    	struct List {
    	    List() : data_(), stride_(0) {}
    	    List(double *data, size_t stride)
		: data_(data), stride_(stride) {}
	    BOOST_GPU_ENABLED
	    double* operator[](size_t i) { return data_ + i*stride_; }
    	private:
    	    double *data_;
	    size_t stride_;
    	};
    	Eri(const Shell::Quartet &quartet,
	    const Context &context,
    	    const rysq::Transpose &transpose);
    	~Eri();
    	void operator()(const Centers &centers,
    			const Quartets &quartets,
    			List eri,
    			const Parameters &parameters,
    			cudaStream_t stream);
	const Shell::Quartet& quartet() const { return quartet_; }
    private:
    	Shell::Quartet quartet_;
    	void *impl_;
    };


} // namespace detail
}
}

#endif /* RYSQ_CUDA_DETAIL_HPP */
