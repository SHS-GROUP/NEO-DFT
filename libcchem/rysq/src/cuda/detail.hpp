#ifndef RYSQ_CUDA_DETAIL_HPP
#define RYSQ_CUDA_DETAIL_HPP

#include <cuda.h>
#include <cuda_runtime.h>

#include <memory>

#include <boost/config.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/bool.hpp>
// #include  <boost/type_traits/is_const.hpp>
#include <boost/array.hpp>
#include <boost/utility/enable_if.hpp>

#include <cmath>

#include <rysq/eri.hpp>
#include <rysq/fock.hpp>
#include "rysq/cuda-memory.hpp"
#include "transpose.hpp"

#include "boost/cuda/runtime.hpp"
#include "boost/cuda/array.hpp"

#include "externals/cxx/array.hpp"
#include "externals/cxx/utility/permute.hpp"
#include "externals/cxx/math/math.hpp"


namespace rysq {
namespace detail {

    struct fock_index {
	template<int I, int J, typename T = int> struct tuple;
	template<int> struct ij;
	template<int> struct kl;
	template<int,int> struct index1;
    };

    template<> struct fock_index::index1<0,1> { static const int value = 0; };
    template<> struct fock_index::index1<2,3> { static const int value = 1; };
    template<> struct fock_index::index1<0,2> { static const int value = 2; };
    template<> struct fock_index::index1<0,3> { static const int value = 3; };
    template<> struct fock_index::index1<1,2> { static const int value = 4; };
    template<> struct fock_index::index1<1,3> { static const int value = 5; };

    template<int I, int J, typename T>
    struct fock_index::tuple {
	static const int value = fock_index::index1<I,J>::value;
	T i, j;
	BOOST_GPU_ENABLED
	tuple() : i(), j() {}
	template<typename U>
	BOOST_GPU_ENABLED
	tuple(const U (&elems)[4]) {
	    this->i = elems[I];
	    this->j = elems[J];
	}
    };
	
    template<> struct fock_index::ij<0> { typedef tuple<0,1> type; };
    template<> struct fock_index::ij<1> { typedef tuple<2,3> type; };
    template<> struct fock_index::ij<2> { typedef tuple<0,2> type; };
    template<> struct fock_index::ij<3> { typedef tuple<0,3> type; };
    template<> struct fock_index::ij<4> { typedef tuple<1,2> type; };
    template<> struct fock_index::ij<5> { typedef tuple<1,3> type; };

    template<> struct fock_index::kl<0> { typedef tuple<2,3> type; };
    template<> struct fock_index::kl<1> { typedef tuple<0,1> type; };
    template<> struct fock_index::kl<2> { typedef tuple<1,3> type; };
    template<> struct fock_index::kl<3> { typedef tuple<1,2> type; };
    template<> struct fock_index::kl<4> { typedef tuple<0,3> type; };
    template<> struct fock_index::kl<5> { typedef tuple<0,2> type; };

}
}

namespace rysq {
namespace cuda {
namespace detail {

    namespace mpl = boost::mpl;
    using rysq::detail::fock_index;

    template<typename T>
    struct matrix_set {
	typedef ushort size_array[4];
	size_array num_blocks, block, base_index;
    private:
	T const* density_[6];
	T *fock_[6];

    public:
	template<class M>
	explicit matrix_set(const M &d, M &f){

	    struct set {
		size_t size[4], block[4];
		set(const M &m) {
		    size[0] = m.get(0,1)->size1();
		    size[1] = m.get(0,1)->size2();
		    size[2] = m.get(2,3)->size1();
		    size[3] = m.get(2,3)->size2();
	    
		    block[0] = m.get(0,1)->block().size1();
		    block[1] = m.get(0,1)->block().size2();
		    block[2] = m.get(2,3)->block().size1();
		    block[3] = m.get(2,3)->block().size2();
		}
	    } d_(d), f_(f);

	    for (int i = 0; i < 6; ++i) {
		density_[i] = d[i]->origin();
		fock_[i] = f[i]->origin();
	    }

	    for (int i = 0; i < 4; ++i) {
		BOOST_VERIFY(f_.size[i] == d_.size[i]);
		BOOST_VERIFY(f_.block[i] == d_.block[i]);
		BOOST_VERIFY(f_.size[i]%f_.block[i] == 0);
		this->num_blocks[i] = f_.size[i]/f_.block[i];
		this->block[i] = f_.block[i];
		BOOST_VERIFY(f.base_index[i] == d.base_index[i]);
		this->base_index[i] = f.base_index[i];
	    }
	}

	T* fock(int i) { return fock_[i]; }
	const T* density(int i) { return density_[i]; }

	template<int I, int J>
	__host__ __device__
	T* fock(const fock_index::tuple<I,J> &index) {
	    return (fock_[fock_index::tuple<I,J>::value] +
		    this->index(index)*size(index));
	}

	template<int I, int J>
	__host__ __device__
	const T* density(const fock_index::tuple<I,J> &index) const {
	    return (density_[fock_index::tuple<I,J>::value] +
		    this->index(index)*size(index));
	}

	template<int I, int J>
	__host__ __device__
	int size(const fock_index::tuple<I,J>&) const {
	    return int(block[I])*int(block[J]);
	}

	template<int I, int J>
	__host__ __device__
	int index(const fock_index::tuple<I,J> &index) const {
	    return ((index.i - base_index[I]) +
		    (index.j - base_index[J])*num_blocks[I]);
	}
    };

    namespace host = boost::cuda;

    struct Quartet : boost::noncopyable {
	Quartet() : initialized_(false) {}
	Quartet(const rysq::Quartet<rysq::Shell> &quartet)
	    : host_(new rysq::Quartet<rysq::Shell>(quartet)),
	      initialized_(false) {}
	operator rysq::Quartet<rysq::type>() const { return *host_; }
	const rysq::Shell& operator[](size_t i) const {
	    return (*host_)[i];
	} 
	size_t size() const { return host_->size(); }
	size_t K() const { return host_->K(); }
	size_t L() const { return host_->L(); }
	host::array<const void> data() const {
	    initialize();
	    return host::array<const void>(data_.begin(), data_.size());
	}
	host::array<const ushort3> index2d() const {
	    initialize();
	    return index2d_;
	}
	void reset() {
	    data_.reset();
	    index2d_.reset();
	    initialized_ = false;
	}
    private:
	std::auto_ptr<rysq::Quartet<rysq::Shell> > host_;
	mutable host::vector<void> data_;
	mutable host::array<const ushort3> index2d_;
	mutable bool initialized_;
	void initialize() const;
    };

    // inline std::ostream& operator<<(std::ostream& stream, const Quartet &q) {
    // 	return boost::operator<<(stream, q);
    // }

    struct Quartets {
	typedef const Int4* const_pointer;
	const Int4 *data_;
	size_t size_;
	Quartets() {}
	Quartets(boost::cuda::vector<Int4> &data)
	    : data_(data.begin().get()), size_(data.size()) {}
	Quartets(device_ptr<const Int4> data, size_t size)
	    : data_(data.get()), size_(size) {}
	const_pointer data()const { return data_; }
	__host__ __device__
	size_t size() const { return size_; }
	__host__ __device__ 
	const Int4& operator[](int index) const {
	    return data_[index];
	}
    };

    struct Centers {
	typedef const rysq::Center* pointer;
	Centers() : data_() {}
	Centers(device_ptr<const rysq::Center> data)
	    : data_(data.get()) {}
	Centers(const boost::cuda::vector<const rysq::Center> &data)
	    : data_(data.begin().get()) {}
	BOOST_GPU_ENABLED
	operator pointer() const { return data_; }
    private:
	pointer data_;
    };

    struct Eri {
	typedef rysq::Eri::Parameters Parameters;
	struct List {
	    List() : data_() {}
	    List(boost::cuda::device_ptr<double*> data)
		: data_(data.get()) {}
	    operator double* const*() { return data_; }
	private:
	    double* const *data_;
	};
	Eri(const rysq::Quartet<Shell> &quartet,
	    const rysq::Transpose &transpose);
	~Eri();
	void operator()(const Centers &centers,
			const Quartets &quartets,
			List eri,
			const Parameters &parameters,
			const boost::cuda::stream &stream);
	void reset() { quartet_.reset(); }
    private:
	Quartet quartet_;
	void *impl_;
    };


    struct Fock {
	//typedef detail::fock_set<double> fock_set;
	// typedef detail:: matrix_set<double, true>  fock_set;
	typedef rysq::Fock::Parameters Parameters;
	typedef detail::matrix_set<double> Set;
	Fock(const rysq::Quartet<Shell> &quartet,
	     const rysq::Transpose &transpose);
	~Fock();
	void operator()(const detail::Centers &centers,
			const detail::Quartets &quartets,
			const Set &set, int* const (&mutex)[6],
			const Parameters &p,
			const boost::cuda::stream &stream);
	void reset() { quartet_.reset(); }
    private:
	Quartet quartet_;
	void *impl_;
    };

} // namespace detail
}
}

#endif /* RYSQ_CUDA_DETAIL_HPP */
