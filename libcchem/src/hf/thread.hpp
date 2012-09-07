#ifndef HF_THREAD_HPP
#define HF_THREAD_HPP

#include <vector>
#include <boost/array.hpp>
// #include <boost/thread.hpp>
// #include <boost/thread/exceptions.hpp>
// #include <boost/thread/mutex.hpp>

#include <boost/typeof/typeof.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include "hf/hf.hpp"
#include "hf/work_queue.hpp"

#include <rysq.hpp>
#include "adapter/rysq.hpp"
#include "utility/timer.hpp"

#include "omp.hpp"

namespace hf {

    using utility::timer;


    struct Transpose {
	static const int BRA = 1;
	static const int KET = 2;
	static const int BRAKET = 4;
	explicit Transpose(int value = 0) {
	    int q[4] = { 0, 1, 2, 3};
	    if (value & BRA) std::swap(q[0], q[1]);
	    if (value & KET) std::swap(q[2], q[3]);
	    if (value & BRAKET) {
		std::swap(q[0], q[2]);
		std::swap(q[1], q[3]);
	    }
	    index_ = 0;
	    for (int i = 0; i < 4; ++i) {
		index_ = (index_ | q[i] << i*2); 
	    }
	}
	template<class Q>
	Q operator()(const Q &q) const {
	    Q p;
	    for (int i = 0; i < 4; ++i) {
		p[i] = q[index_ >> i*2 & 0x3];
	    }
	    return p;
	}
    private:
	int index_;
    };


    template<typename T>
    class lockable : boost::noncopyable {
    public:
	typedef T& reference;
	typedef const T& const_reference;
	explicit lockable(T &data, const boost::array<size_t,2> &shape)
	    : data_(data),
	      mutex_(new omp::lock[shape[0]*shape[1]]),
	      array_(mutex_, shape) {}
	~lockable() { delete[] mutex_; }
	omp::lock& lock(int i, int j) const {
	    return array_[j][i];
	}
	reference data(){return data_; }
	const_reference data()const {return data_; }
    private:
	T &data_;
	mutable omp::lock *mutex_;
	mutable boost::multi_array_ref<omp::lock,2> array_;
    };

    template<typename T>
    lockable<T> make_lockable(T &data) {
	return lockable <T>(data);
    }

    template<class Matrix, class Quartets = std::vector<boost::array<int,4> > >
    class Thread : boost::noncopyable {
    public:

	struct Device;

	struct {
	    Matrix matrix;
	    int index[2];
	} fock_[6];

	template<class Blocks, class Centers, class Density, class Fock,
		 class Screen>
	void operator()(const Blocks &blocks, const Centers &centers,
			const Density &D, lockable<Fock> &F,
			const Screen &screen, work_queue &queue,
			Parameters parameters) {
	    timer t;
	    work_queue::value_type task;
	    Quartets quartets;
	    while (true) {
		try { task = queue.pop(); }
		catch (const std::exception &e) { break; }
		std::swap(task[1], task[2]); // ikjl -> ijkl
		BOOST_AUTO(S, (screen(blocks, task, 1<<16)));
		while (true) {
		    S(quartets);
		    if (quartets.empty()) break;
		    this->operator()(blocks, centers, D, F,
				     task, quartets, parameters);
		    quartets.clear();
		}
	    }
	    // std::cout << "CPU thread time: " << t << std::endl;
	}


	template<class Blocks, class Centers, class Density, class Fock,
		 class Task>
	void operator()(const Blocks &blocks, const Centers &centers,
			const Density &D, lockable<Fock> &F,
			const Task &task, const Quartets &quartets,
			const Parameters &parameters) {

	    if (quartets.empty()) return;

	    const Basis::Block* block[4];
	    int base[4];

	    for (int i = 0; i < 4; ++i) {
		block[i] = &blocks[task[i]];
		base[i] = block[i]->start();
	    }

	    rysq::hf::matrix_set<adapter::rysq::density_matrix> density(base);
	    rysq::hf::matrix_set<adapter::rysq::fock_matrix> fock(base);

	    typedef adapter::rysq::Shell Shell;
	    typedef rysq::Quartet<rysq::Shell> Quartet;
	    Shell a(block[0]->shell());
	    Shell b(block[1]->shell());
	    Shell c(block[2]->shell());
	    Shell d(block[3]->shell());
	    Quartet quartet(a,b,c,d);
	    
	    int k = 0;
	    foreach (const size_t (&index)[2],  rysq::index_list()) {
	        size_t i = index[0], j = index[1];
	        size_t im = std::distance(&blocks[0], block[i]);
	        size_t jm = std::distance(&blocks[0], block[j]);

	        density.get(i,j) = adapter::rysq::density_matrix(D.m(im,jm));
	        density.get(j,i) = adapter::rysq::density_matrix(D.m(jm,im));

		{
		    typename Fock::matrix &f = F.data().m(im,jm);
		    fock_[k].matrix.resize(f.size1(), f.size2(),
					   f.block1(), f.block2());
		    fock_[k].matrix.clear();
		    fock_[k].index[0] = im;
		    fock_[k].index[1] = jm;
		    fock.get(i,j) = adapter::rysq::fock_matrix(fock_[k].matrix);
		    ++k;
		}
	    }

	    rysq::Fock f(quartet);
	    f(centers, quartets, density, fock,
	      rysq::Fock::Parameters(parameters.cscale,
				     parameters.xscale,
				     parameters.cutoff));

	    for (int k = 0; k < 6; ++k) {
		int i = fock_[k].index[0], j = fock_[k].index[1];
		omp::scoped_lock lock(F.lock(i,j));
		F.data().m(i,j) += fock_[k].matrix;
	    }

	}
    };

}

#endif /* HF_THREAD_HPP */
