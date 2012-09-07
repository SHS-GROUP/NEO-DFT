#ifndef TRANSFORM_HPP
#define TRANSFORM_HPP

#include "core/integral.hpp"
#include "blas.hpp"

#include <memory>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/mpl/if.hpp>
#include <boost/thread/thread.hpp>
#include "boost/threadpool.hpp"
#include <boost/shared_ptr.hpp>

#include <boost/ref.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/ptr_container/ptr_vector.hpp>


#include "tensor/array.hpp"
#include "externals/boost/parallel.hpp"
#include "boost/numeric/ublas/storage_adaptors.hpp"

#ifdef HAVE_GPU_INTEGRALS
#define HAVE_GPU_TRANSFORM
#include "core/integral/gpu.hpp"
#include "boost/cuda/runtime.hpp"
#include "boost/cuda/exception.hpp"
#include "boost/numeric/bindings/cublas/cublas.hpp"
#include "gpu/multi_array.hpp"
#include "multi_array/gpu/algorithm.hpp"
#endif

#include "runtime.hpp"

namespace  transform {

    typedef boost::numeric::ublas::column_major layout;
    typedef boost::numeric::ublas::upper upper;
    typedef boost::numeric::ublas::reserve_array<double> storage;
    typedef boost::numeric::ublas::matrix<double, layout, storage> Matrix;
    typedef boost::numeric::ublas::array_adaptor<double> array_type;
    typedef boost::numeric::ublas::matrix<
	double, layout, array_type> matrix_adapter;
    typedef boost::numeric::ublas::triangular_matrix<
	double, upper, layout, array_type> triangular_matrix_adapter;
    typedef Integral::Order Order;    

    struct apply;

    template<class Transform>
    struct Lambda;

#ifndef HAVE_GPU_TRANSFORM
    struct Gpu {
	explicit Gpu(int) {}
	static void reset() {}
    };
#else
    struct Gpu : boost::cuda::thread {
	// typedef boost::cuda::thread::disabled disabled;
	// struct Device {
	//     Device() : enabled_(false) {}
	//     explicit Device(int device) : enabled_(true), device_(device) {}
	//     int get() const { return device_; }
	//     operator bool() const { return enabled_; }
	// private:
	//     bool enabled_;
	//     int device_;
	// };
	typedef boost::numeric::bindings::cublas::matrix<double> Matrix;
	typedef gpu::multi_array<double,4> Array1;
	typedef Integral::Gpu Eri;
	explicit Gpu(int device) :
	    boost::cuda::thread(device, boost::cuda::flags::map_host) {}
	static void reset() { boost::cuda::reset(); }
	Matrix& C1() { return C1_; }
	Matrix& matrix1() { return m1_; }
	Matrix& matrix2() { return m2_; }
	Matrix& matrix3() { return m3_; }
	Array1& array1() { return array1_; }
	operator bool() const { return boost::cuda::thread::enabled(); }
	Eri& eri() { return eri_; }
    private:
	Integral::Gpu integrals_;
	Matrix C1_, m1_, m2_, m3_;
	Array1 array1_;
	Eri eri_;
    public:
	void synchronize() { boost::cuda::thread::synchronize(); }
	typedef boost::array<Integral::index,4> quartet;
	boost::cuda::vector<quartet> quartets;
	boost::cuda::vector<double*> data;
	multi_array::gpu::algorithm::pack pack;
    };
#endif

    template<bool triangular>
    struct Transform : boost::noncopyable {

	typedef typename boost::mpl::if_c<triangular,
					  Transform,
					  Transform<true> >::type Triangular;

	static const size_t N = (triangular) ? 3 : 4;
	typedef tensor_array<4,double> Array1;
	typedef tensor_array<N,double> Array2;

	typedef transform::Matrix Matrix;
	typedef typename boost::mpl::if_c<triangular,
					  transform::triangular_matrix_adapter,
					  transform::matrix_adapter>::type Matrix2;
	typedef transform::Order Order;    
	struct Evaluate;

	//Lambda<Transform> operator()();

    private:
	struct Slave {
	    explicit Slave(int id) : id(id) {}
	    int id;
	};
	Transform(const Transform &master, Slave device)
	    : order_(master.order_),
	      gpu_(), device_(device.id) {}

    public:
	Transform(Order order = Order(), runtime rt = runtime())
	    : order_(order), gpu_()
	{
	    size_t num_threads = rt.get<int>("threads");
	    std::vector<int> devices = rt.get< std::vector<int> >("gpu::devices");

#ifndef HAVE_GPU_TRANSFORM
	    devices.clear();
#endif
	    num_threads = std::max(devices.size(), num_threads) - devices.size();
	    for (size_t i = 1; i < num_threads; ++i) {
	        threads_.push_back(new Transform(*this, Slave(-1)));
	    }

	    if (rt.profile()) {
		std::cout <<  "cpu threads: " << num_threads;
		if (!devices.empty()) {
		    std::cout << ", gpu devices: ";
		    foreach (int i, devices) {
			std::cout << i << " ";
		    }
		}
		std::cout << std::endl;
	    }

	    if (!devices.empty()) {
	    	// first gpu on this thread
	    	Gpu::reset();
	    	device_ = devices.back();
	    	devices.pop_back();
	    	gpu_.reset(new Gpu(device_));
	    }
	    // foreach (int i, devices) {
	    // 	//threads_.push_back(new Transform(*this, Slave(i)));
	    // }


	    pool_ = boost::threadpool::pool(threads_.size());
	}

	const Order& order() const { return order_; }

	template<typename U>
	void resize(const U (&dims)[4]) {
	    boost::array<size_t,N> dims_;
	    if (N == 3) {
		if (dims[0] != dims[1])
		    throw std::range_error("invalid triangular dimensions");
		size_t n = dims[0];
		dims_[0] = (n*n + n)/2;
		std::copy(dims + 2, dims + 4, &dims_[1]);
	    }
	    else {
		std::copy(dims, dims + 4, &dims_[0]);
	    }
	    array2_.resize(dims_);
	}

	Array1& array1() { return array1_; }

	Array2& array2() { return array2_; } 
	const Array2& array2() const { return array2_; } 

	operator double*() { return array2_.data(); }
	operator const double*() const { return array2_.data(); }

	Matrix& matrix1() { return matrix_[0]; }
	Matrix& matrix2() { return matrix_[1]; }
	Matrix& matrix3() { return matrix_[2]; }	    
	
	template<class F>
	void operator()(const Matrix &C1, const Matrix &C2,
			integral::generator<1> integrals, F f) {
	    blas::set_num_threads(1);

	    boost::parallel parallel;
	    boost::thread_group threads;
	    BOOST_AUTO(shells, parallel.range(integrals.basis().shells()));

	    foreach (Transform &T, threads_) {
	     	using boost::ref;
		using boost::bind;
		pool_.schedule(bind<void>(ref(T),
		 			  ref(C1), ref(C2),
					  ref(integrals), f, ref(shells)));
		// threads.add_thread(new boost::thread(ref(T),
		// 				     ref(C1), ref(C2),
		// 				     ref(integrals),
		// 				     f, ref(shells)));
	    }
	    apply(C1, C2, integrals, f, shells);
	    pool_.wait(0);
	    threads.join_all();
	}

	template<class S, class F>
	void operator()(const Matrix &C1, const Matrix &C2,
			const integral::generator<1> &integrals,
			const F &f, const S &shells) {
	    if (device_ > -1) gpu_.reset(new Gpu(device_));
	    apply(C1, C2, integrals, f, shells);
	    gpu_.reset();
	}

	struct Functor {
	    template<int> struct index {};
	    Functor(Transform &T) : T(T) {}
	    template<class R, class E, class F>
	    void operator()(index<4>, const R &range, const Array1 &A,
			    const boost::numeric::ublas::matrix_expression<E> &C,
			    const F &f) const {
		BOOST_AUTO(begin, range.data().front());
		foreach (int k, range) {
		    int j = int((1 + sqrt(8*k +1))/2 - 1);
		    int i = k - (j*j + j)/2;
		    // std::cout << "xxx: " << k << " " << i << " " << j << std::endl;
		    BOOST_AUTO(&a, T.matrix1());
		    BOOST_AUTO(&t, T.matrix2());
		    size_t start[] = { 0, k-begin, 0, 0 };
		    size_t stop[] = { A.size<0>(), k-begin+1, A.size<2>(), 1 };
		    a.resize(stop[0], stop[2], false);
		    A.get(a.data().begin(), start, stop);
		    T.gemm(a, C, t);
		    f(i, j, t);
		}
	    }
	    Transform &T;
	};

	template<class E, class F>
	void operator()(size_t begin, size_t end, const Array1 &A,
			const boost::numeric::ublas::matrix_expression<E> &C,
			const F &f) {
	    blas::set_num_threads(1);

	    std::vector<int> indices;
	    for (size_t k = begin; k < end; ++k) {
		indices.push_back(k);
	    }

	    boost::parallel parallel;
	    BOOST_AUTO(r, parallel.range(indices));

	    typename Functor::template index<4> index;
	    foreach (Transform &T, threads_) {
	    	using boost::bind;
		using boost::ref;
	    	pool_.schedule(bind<void>(Functor(T), index,
	    				  ref(r), ref(A), ref(C), f));
	    }
	    (Functor(*this))(index, r, A, C, f);
	    pool_.wait(0);
	}
	
	template<class E0, class E1, class F>
	typename boost::disable_if<boost::is_same<Matrix,F> >::type
	operator()(const boost::numeric::ublas::matrix_expression<E0> &a,
		   const boost::numeric::ublas::matrix_expression<E1> &b,
		   const F &f) {
	    Matrix &T2 = matrix1();
	    blas::set_num_threads(1+threads_.size());
	    gemm(a, b, T2);
	    f(T2);
	}

	void operator()(const Matrix &C1, const Matrix &C2,
			const integral::generator<2> &integrals);

    public:
	typedef transform::Gpu Gpu;
	Gpu* gpu() { return gpu_.get(); };

    protected:
	Order order_;
	Array1 array1_;
	Array2 array2_;
	Matrix matrix_[3];
	integral::array integrals_;
	boost::ptr_vector<Transform> threads_;
	std::auto_ptr<Gpu> gpu_;
	int device_;
	boost::mutex mutex_;
	boost::threadpool::pool pool_;

    private:
	template<class S, class F>
	void apply(const Matrix &C1, const Matrix &C2,
		   const integral::generator<1> &integrals,
		   const F &f, const S &shells) {
	    BOOST_PROFILE_FUNCTION();
#ifdef HAVE_GPU_TRANSFORM
	    if (gpu()) {
		gpu()->eri().set(integrals.basis().centers());
		gpu()->C1() = C1;
	    }
#endif
	    typedef integral::generator<1>::Shell Shell;
	    foreach (const Shell &P, shells) {
		operator()(C1, C2, integrals.push_front(P));
		f(*this, P);
	    }
	}

	template<class A, class B>
	void gemm(const boost::numeric::ublas::matrix_expression<A> &a,
		  const boost::numeric::ublas::matrix_expression<B> &b,
		  Matrix &c, size_t num_threads = 1) {
	    c.resize(a().size1(), b().size2(), false);
	    blas::gemm(1, a(), b(), 0, c);
	}

    };


}

#include "core/transform/contract2.hpp"

namespace transform {
    
    template<bool Symmetry>
    void Transform<Symmetry>::operator()(const Matrix &C1, const Matrix &C2,
					 const integral::generator<2> &integrals) {
	size_t size[] = { C1.size1(), C2.size1(),
			  integrals.get<0>().size(), integrals.get<1>().size() };
	this->resize(size);
	contract2(C1, C2, integrals, *this);
    }

}

typedef transform::Transform<false> Transform;


#endif // TRANSFORM_HPP
