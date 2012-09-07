#ifndef BOOST_NUMERIC_BINDINGS_CUBLAS_STORAGE_HPP
#define BOOST_NUMERIC_BINDINGS_CUBLAS_STORAGE_HPP

#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/bindings/ublas/storage.hpp>
//#include <boost/numeric/bindings/traits/vector_traits.hpp>
#include <cuda_runtime.h>
#include <cublas.h>

#include <boost/numeric/bindings/cublas/exception.hpp>


namespace boost {
namespace numeric {
namespace bindings {
namespace cublas {

    template<typename T>
    struct array_adaptor :
	ublas::storage_array<array_adaptor<T> >,
	private ublas::array_adaptor<T>
    {
	typedef array_adaptor<T> self_type;
	// friend class traits::default_vector_traits<self_type, T>;
	// friend class traits::default_vector_traits<const self_type, T>;

	typedef ublas::array_adaptor<T> ublas_base;
	using ublas_base::size_type;
	using ublas_base::difference_type;
	typedef typename ublas_base::pointer pointer;
	//using ublas_base::pointer;
	using ublas_base::const_pointer;
	using ublas_base::iterator;
	using ublas_base::const_iterator;
	using ublas_base::value_type;

	array_adaptor() : ublas_base(0, pointer(0)) {}
	array_adaptor(size_t size, pointer data) : ublas_base(size, data) {}
	void swap(array_adaptor &a) {
	    ublas_base::swap(a);
	}

    public:
	using ublas_base::size;
	using ublas_base::begin;
	using ublas_base::end;
	using ublas_base::rbegin;
	using ublas_base::rend;
    };
    
    template<typename T>
    struct unbounded_array :
	ublas::storage_array<unbounded_array<T> >,
	private ublas::array_adaptor<T>
    {
	typedef unbounded_array<T> self_type;
	// friend class traits::default_vector_traits<self_type, T>;
	// friend class traits::default_vector_traits<const self_type, T>;

	typedef ublas::array_adaptor<T> ublas_base;
	using ublas_base::size_type;
	using ublas_base::difference_type;
	using ublas_base::pointer;
	using ublas_base::const_pointer;
	using ublas_base::iterator;
	using ublas_base::const_iterator;
	using ublas_base::value_type;
	unbounded_array(size_t size = 0)
	    : ublas_base(size, allocate(size)) {}
	unbounded_array(const unbounded_array &data)
	    : ublas_base(0, (T*)0)
	{
	    assign(data);
	}
	unbounded_array& operator=(const unbounded_array &data) {
	    assign(data);
	}
	~unbounded_array() { free(this->begin()); }
	void resize(size_t size) {
	    if (size == ublas_base::size()) return;
	    {
		unbounded_array null(0);
		this->swap(null);
	    }
	    unbounded_array data(size);
	    this->swap(data);
	}
	void swap(unbounded_array &a) {
	    ublas_base::swap(a);
	}
    public:
	using ublas_base::size;
	using ublas_base::begin;
	using ublas_base::end;
	using ublas_base::rbegin;
	using ublas_base::rend;
    private:
	static void throw_(cudaError_t error) {
	    if (error != cudaSuccess) {
		throw std::runtime_error(cudaGetErrorString(error));
	    }
	}
	void assign(const unbounded_array &data) {
	    resize(data.size());
	    if (!data.size()) return;
	    throw_(cudaMemcpy(this->begin(), data.begin(), data.size()*sizeof(T),
			      cudaMemcpyDeviceToDevice));
	}
	static void free(T *data){
	    if (data) {
		throw_(cudaFree(data));
		// cublasFree(data);
		// check_status();
	    }
	}
	static T* allocate(size_t size) {
	    void* data = NULL;
	    if (size) {
		throw_(cudaMalloc(&data, size*sizeof(T)));		
		// cublasAlloc(size, sizeof(T), &data);
		// check_status();
	    }
	    return static_cast<T*>(data);
	}
    };

}
}
}
}

#endif // BOOST_NUMERIC_BINDINGS_CUBLAS_STORAGE_HPP
