#ifndef RYSQ_CUDA_HPP
#define RYSQ_CUDA_HPP


#include <vector>
#include <memory>
#include <boost/noncopyable.hpp>

#include "rysq-core.hpp"
#include "rysq-eri.hpp"
#include "rysq-fock.hpp"

extern "C" {
    struct cudaDeviceProp;
}

namespace rysq {
namespace cuda {


    namespace runtime {

	struct exception: std::runtime_error {
	    explicit exception(const char *what)
		: std::runtime_error(what) {}
	};

	struct Stream {
	    explicit Stream(void *data = NULL) : data_(data) {}
	    void synchronize();
	    void* data() const { return data_; }
	private:
	    void *data_;
	};

	struct Index;

	struct Context {
	    Context();
	    ~Context();
	    const cudaDeviceProp& properties() const;
	    const void* index(rysq::type A, rysq::type B,
			      rysq::type C, rysq::type D) const;
	private:
	    mutable std::auto_ptr<Index> index_;
	    mutable std::auto_ptr<cudaDeviceProp> properties_;
	};


	void* malloc(size_t size);
	void free(void *ptr);

	void memset(void *ptr, int value, size_t size);


	typedef enum {
	    device_to_device = 1,
	    device_to_host = 2,
	    host_to_device = 4 
	} copy_type;

	template<typename T>
	void copy(size_t size, const T *from, T *to, copy_type t);

	template<>
	void copy(size_t size, const void *from, void *to, copy_type);

	template<typename T>
	void copy(size_t size, const T *from, T *to, copy_type t) {
	    copy(size*sizeof(T), (const void*)from, (void*)to, t);
	}


	template<typename T>
	void copy(size_t size, const T *from, T *to, copy_type t,
		  const Stream &s);

	template<>
	void copy(size_t size, const void *from, void *to, copy_type,
		  const Stream &s);

	template<typename T>
	void copy(size_t size, const T *from, T *to, copy_type t,
		  const Stream &s) {
	    copy(size*sizeof(T), (const void*)from, (void*)to, t, s);
	}


    }


    namespace detail {

	template<typename T>
	struct device_vector {
	    typedef T value_type;

	    device_vector(size_t size = 0)
		: data_(), size_(), capacity_()
	    {
		resize(size);
	    }

	    ~device_vector() {
		runtime::free(data_);
	    }

	    device_vector(const device_vector &o)
		: data_(), size_(), capacity_()
	    {
		assign(o.size(), o.begin(), runtime::device_to_device);
	    }

	    void operator=(const device_vector &o) {
		assign(o.size(), o.begin(), runtime::device_to_device);
	    }

	    void resize(size_t size) {
		if (size > capacity_) {
		    runtime::free(data_);
		    data_ = static_cast<T*>(runtime::malloc(size*sizeof(T)));
		    capacity_ = size;
		}
		size_ = size;
	    }

	    T* begin() { return data_; }
	    const T* begin() const { return data_; }
	    size_t size() const { return size_; }

	    void assign(size_t size, const T *data, runtime::copy_type t) {
		resize(size);
		if (!size) return;
		copy(size, data, data_, t);
	    }

	    void assign(size_t size, const T *data, runtime::copy_type t,
			const runtime::Stream &s) {
		resize(size);
		if (!size) return;
		copy(size, data, data_, t, s);
	    }

	private:
	    T *data_;
	    size_t size_, capacity_;
	};

    }


    struct Centers : detail::device_vector<Center> {
	typedef detail::device_vector<Center> vector;
	void assign(const std::vector<Center> &v) {
	    vector::assign(v.size(), &v[0], runtime::host_to_device);
	}
    };

    struct Quartets : detail::device_vector<Int4> {
	typedef detail::device_vector<Int4> vector;
	void assign(const std::vector<Int4> &v, const runtime::Stream &s) {
	    vector::assign(v.size(), &v[0], runtime::host_to_device, s);
	}
    };

    struct Mutex : private detail::device_vector<int> {
	typedef detail::device_vector<int> vector;
	void resize(size_t size) {
	    size_ = size;
	    vector::resize(size*size);
	    runtime::memset(this->data(), 0, size*size*sizeof(int));
	}
	int* data() { return vector::begin(); }
	size_t size() const { return size_; }
    private:
	size_t size_;
    };


    template<typename T, class A = detail::device_vector<T> >
    struct block_matrix : block_matrix_base<T> {
	typedef A array_type;
	typedef block_matrix_base<T> base;
	using base::size1;
	using base::size2;
	using base::size;
	using base::block;
	block_matrix() {}
	T* block(int a,int b) {
	    return data_.begin() +  base::layout_.block_at(a,b);
	}
	const T* block(int a,int b) const {
	    return data_.begin() +  base::layout_.block_at(a,b);
	}
	void clear() {
	    runtime::memset(data_.begin(), 0, size()*sizeof(T));
	}
	template<class M>
	void assign(const M &a) {
	    base::check_size(a);
	    rysq::block_matrix<T>
		h(size1(), size2(), block().size1(), block().size2());
	    h.assign(a);
	    data_.assign(h.size(), h.data(), runtime::host_to_device);
	}
	void resize(int size1, int size2, int block1, int block2) {
	    base::layout_ = block_matrix_layout(size1, size2, block1, block2);
	    base::size1_ = size1;  base::size2_ = size2;
	    data_.resize(size());
	}
	array_type& data() { return data_; }
	const array_type& data() const { return data_; }
    private:
	array_type data_;
    };


    template<typename T, class A>
    rysq::block_matrix<T>
    host(const cuda::block_matrix<T,A> &m) {
	rysq::block_matrix<T>
	    h(m.size1(), m.size2(), m.block().size1(), m.block().size2());
	namespace rt = runtime;
	rt::copy(h.size(), m.data().begin(), h.data(), rt::device_to_host);
	return h;
    }


    namespace detail {

	class Eri;
	class Fock;

	struct Shell {

	    struct Quartet;
	    struct Primitive {
		double a, C, Cs, _;
	    };

	    BOOST_GPU_ENABLED
	    Shell() : data_(), type_(), K_() {}

	    BOOST_GPU_ENABLED
	    Shell(const Primitive *data, rysq::type type, int K)
		: data_(data), type_(type), K_(K) {}

	    BOOST_GPU_ENABLED
	    const Primitive* data() const { return data_; }

	    BOOST_GPU_ENABLED
	    int K() const { return K_; }

	    BOOST_GPU_ENABLED
	    int L() const { return abs(type_); }

	    BOOST_GPU_ENABLED
	    int nc() const { return (1 + (type_ < 0)*L()); } 

	    BOOST_GPU_ENABLED
	    rysq::type type() const { return rysq::type(type_); }

	    BOOST_GPU_ENABLED
	    operator rysq::type() const { return this->type(); }

	    BOOST_GPU_ENABLED
	    int sp() const { return (this->type() == rysq::SP); }

	    BOOST_GPU_ENABLED
	    int begin() const { return rysq::shell::begin(*this); }

	    // BOOST_GPU_ENABLED
	    // const double& operator()(int i, int j) const {
	    // 	return data_[i + (j+1)*K_];
	    // }

	    BOOST_GPU_ENABLED
	    const Primitive& operator()(int i) const {
		return data_[i];
	    }

	    BOOST_GPU_ENABLED
	    int size() const {
		rysq::type type = *this;
		return (rysq::shell::end(type) - rysq::shell::begin(type));
	    }

	protected:
	    const Primitive *data_;
	    int type_;
	    int K_;
	};

	struct Shell::Quartet {

	    Shell elems[4];
	    char order[4];

	    struct Swap {
		bool bra, ket, braket;
	    } swap_;

	    Quartet(const Shell &a,const Shell &b,
		    const Shell &c, const Shell &d) {
		elems[0] = a;
		elems[1] = b;
		elems[2] = c;
		elems[3] = d;
		for (int i = 0; i < 4; ++i) {
		    order[i] = i;
		}
		swap_.bra = 0;
		swap_.ket = 0;
		swap_.braket = 0;
	    }

	    BOOST_GPU_ENABLED
	    Shell& operator[](int i) { return elems[i]; }

	    BOOST_GPU_ENABLED
	    const Shell& operator[](int i) const { return elems[i]; }

	    BOOST_GPU_ENABLED
	    int size() const {
		return (elems[0].size()*elems[1].size()*
			elems[2].size()*elems[3].size());
	    }

	    BOOST_GPU_ENABLED
	    int L() const {
		return (elems[0].L() + elems[1].L() +
			elems[2].L() + elems[3].L());
	    }

	    BOOST_GPU_ENABLED
	    int K() const {
		return (elems[0].K()*elems[1].K()*
			elems[2].K()*elems[3].K());
	    }

	    BOOST_GPU_ENABLED
	    int sp() const {
		return (elems[0].sp() + elems[1].sp() +
			elems[2].sp() + elems[3].sp());
	    }

	    template<class S>
	    void swap(const S &s) {
		if (s.bra) {
		    swap_.bra = s.bra;
		    swap(0,1);
		}
		if (s.ket) {
		    swap_.ket = s.ket;
		    swap(2,3);
		}
		if (s.braket) {
		    swap_.braket = s.braket;
		    swap(0,2);
		    swap(1,3);
		}
	    }

	    void sort() {
		// permute s.t. S shells are right-most
		if (!elems[0].L())
		    swap(0,1);
		if (!elems[2].L())
		    swap(2,3);
		if (!(elems[0].L() || elems[1].L())) {
		    swap(0,2);
		    swap(1,3);
		}
	    }

	protected:
	    void swap(int i, int j) {
		std::swap(this->order[i], this->order[j]);
		std::swap(this->elems[i], this->elems[j]);
	    }
	};

    }

    typedef runtime::exception exception;
    typedef runtime::Stream Stream;
    typedef runtime::Context Context;

    struct Shell : detail::Shell, boost::noncopyable {
	Shell(const rysq::Shell &shell);
	~Shell();
    };

    struct Eri : boost::noncopyable {
    	typedef rysq::Eri::Parameters Parameters;
    	Eri(Shell &a, Shell &b, Shell &c, Shell &d,
	    const Context &context);
    	~Eri();
    	operator bool() const { return (kernel_.get() != NULL); }
    	void operator()(const Centers &centers,
    			const Quartets &quartets,
    			double *I,
    			const Parameters &parameters = Parameters(),
			Stream s = Stream());
    private:
	Transpose transpose_;
    	std::auto_ptr<detail::Eri> kernel_;
    };


    typedef block_matrix<double> density_matrix;
    typedef block_matrix<double> fock_matrix;


    struct Fock : boost::noncopyable {
	typedef hf::matrix_ptr_set<fock_matrix> fock_matrix_set;
	typedef hf::matrix_ptr_set<density_matrix> density_matrix_set;
	typedef rysq::Fock::Parameters Parameters;
	Fock(Shell &a, Shell &b, Shell &c, Shell &d,
	     const Context &context);
        ~Fock();
	operator bool() const { return kernel_.get(); }
	void operator()(const Centers &centers,
			const Quartets &quartets,
			density_matrix_set D, fock_matrix_set F,
			Mutex &mutex,
			const Parameters &parameters,
			Stream s = Stream());
    private:
	Transpose transpose_;
	boost::array<size_t,4> block_;
	std::auto_ptr<detail::Fock> kernel_;
    };
	    

}
}

#endif /* RYSQ_CUDA_HPP */
