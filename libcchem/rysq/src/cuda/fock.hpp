#ifndef RYSQ_CUDA_FOCK_HPP
#define RYSQ_CUDA_FOCK_HPP


#include "rysq-eri.hpp"
#include "rysq-fock.hpp"
#include "rysq-cuda.hpp"

#include "transpose.hpp"

#include "cuda/detail.hpp"


#include <memory>
#include <boost/config.hpp>

#include <cuda.h>
#include <cuda_runtime.h>


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
	static const int first = I;
	static const int second = J;
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

    // struct Index {
    // 	static void initialize();
    // 	static void finalize();
    // 	static boost::cuda::array<const ushort3> get();
    // };

    using rysq::detail::fock_index;

    template<typename T>
    struct matrix_set {
	int block[4];
	int base_index[4];
	ushort num_blocks[4];
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
		density_[i] = d[i]->data().begin();
		fock_[i] = f[i]->data().begin();
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

	BOOST_GPU_ENABLED
	T* fock(int i) { return fock_[i]; }

	BOOST_GPU_ENABLED
	const T* density(int i) { return density_[i]; }

	template<int I, int J>
	BOOST_GPU_ENABLED
	T* fock(const fock_index::tuple<I,J> &index) {
	    int ij = this->index(index)*size(index);
	    return fock_[fock_index::tuple<I,J>::value] + ij;
	}

	template<int I, int J>
	BOOST_GPU_ENABLED
	const T* density(const fock_index::tuple<I,J> &index) const {
	    int ij = this->index(index)*size(index);
	    return density_[fock_index::tuple<I,J>::value] + ij;
	}

	template<int I, int J>
	BOOST_GPU_ENABLED
	ushort size(const fock_index::tuple<I,J>&idx) const {
	    return block[I]*block[J];
	}

	template<int I, int J>
	BOOST_GPU_ENABLED
	int index(const fock_index::tuple<I,J> &index) const {
	    int i = index.i - base_index[I];
	    int j = index.j - base_index[J];
	    return (i + j*num_blocks[I]);
	}
    };


    struct Mutex {
	explicit Mutex(cuda::Mutex &mutex)
	    : data_(mutex.data()), size_(mutex.size()) {}
	BOOST_GPU_ENABLED
	int* operator[](int i) {
	    return data_ + i*size_;
	}
    private:
	int* data_;
	int size_;
    };


    struct Fock {
	typedef rysq::Fock::Parameters Parameters;
	typedef detail::matrix_set<double> Set;
	Fock(const Shell::Quartet &quartet,
	     const Context &context,
	     const rysq::Transpose &transpose);
	~Fock();
	void operator()(const detail::Centers &centers,
			const detail::Quartets &quartets,
			const Set &set,
			detail::Mutex mutex,
			const Parameters &p,
			cudaStream_t stream);
    private:
	Shell::Quartet quartet_;
	void *impl_;
    };

} // namespace detail
}
}

#endif /* RYSQ_CUDA_FOCK_HPP */
