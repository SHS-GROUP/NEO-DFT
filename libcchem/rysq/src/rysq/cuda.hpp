#ifndef _RYSQ_CUDA_HPP_
#define _RYSQ_CUDA_HPP_

#include <vector>
#include <memory>
#include <limits>
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
// #include <boost/numeric/ublas/matrix.hpp>
#include "rysq/core.hpp"
#include "rysq/eri.hpp"
#include "rysq/fock.hpp"


#include "rysq/cuda-memory.hpp"
#include "rysq/cuda-matrix.hpp"

#include "boost/cuda/vector.hpp"

namespace rysq {
namespace cuda {

    std::vector<int> get_devices();

    void initialize(int device);
    void finish();

    struct thread {
	thread(int device);
	~thread();
    };

    namespace detail {
	class Eri;
	class Fock;
    }

    typedef boost::cuda::vector<Center> Centers;
    typedef std::vector<Int4, boost::cuda::mapped_allocator<Int4> > Quartets;

    struct Eri : boost::noncopyable {
	typedef std::vector<Int4> Quartets;
	typedef boost::cuda::stream stream_type;
	typedef rysq::Eri::Parameters Parameters;
	Eri(const rysq::Quartet<rysq::Shell> &shells,
	    const stream_type &stream = synchronous);
	~Eri();
	operator bool() const { return (kernel_.get() != NULL); }
	void operator()(const Centers &centers,
			const Quartets &quartets,
			const std::vector<double*> &I,
			const Parameters &parameters = Parameters());
	void operator()(const Centers &centers, size_t size,
			boost::cuda::device_ptr<const Int4> quartets,
			boost::cuda::device_ptr<double*> I,
			const Parameters &parameters = Parameters());
	stream_type& stream() { return stream_; }
    private:
	stream_type stream_;
	std::auto_ptr<detail::Eri> kernel_;
	boost::cuda::vector<Int4> quartets_;
	boost::cuda::vector<double*> integrals_;
    };

    typedef cuda::block_matrix<double> density_matrix;
    typedef cuda::block_matrix<double> fock_matrix;

    template<typename T, class A>
    rysq::block_matrix<T, boost::shared_array<T> >
    host(const block_matrix<T,A> &m) {
	rysq::block_matrix<T,boost::shared_array<T> >
	    h(m.size1(), m.size2(), m.block().size1(), m.block().size2());
	boost::cuda::copy(m.data().begin(), h.data(), h.size());
	return h;
    }


    struct Fock : boost::noncopyable {
	//typedef cuda::Quartets Quartets;
	typedef std::vector<Int4> Quartets;
	typedef hf::matrix_ptr_set<fock_matrix> fock_matrix_set;
	typedef hf::matrix_ptr_set<density_matrix> density_matrix_set;
	typedef rysq::Fock::Parameters Parameters;
	Fock();
	Fock(const rysq::Quartet<rysq::Shell> &shells);
	~Fock();
	void swap(Fock &other);
	operator bool() const { return kernel_.get() != NULL; }
	void operator()(const cuda::Centers &centers,
			const Quartets &quartets,
			density_matrix_set D, fock_matrix_set F,
			const Parameters &parameters);
    private:
	Transpose transpose_;
	std::auto_ptr<detail::Fock> kernel_;
	boost::array<size_t,4> block_;
	boost::shared_ptr<boost::cuda::vector<int> > mutex_;
	boost::cuda::vector<Int4> quartets_;
	std::auto_ptr<boost::cuda::stream> stream_;
    };
	    

}
}

#endif /* _RYSQ_CUDA_HPP_ */
