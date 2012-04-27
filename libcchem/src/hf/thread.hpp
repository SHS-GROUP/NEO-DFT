#ifndef HF_THREAD_HPP
#define HF_THREAD_HPP

#include <vector>
#include <boost/array.hpp>
#include <boost/thread.hpp>
#include <boost/thread/exceptions.hpp>
#include <boost/thread/mutex.hpp>

#include <boost/typeof/typeof.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include "hf/hf.hpp"

#include <rysq.hpp>
#include "adapter/rysq.hpp"
#include "utility/timer.hpp"

namespace hf {

    using utility::timer;


    template<typename T>
    class lockable : boost::noncopyable {
    public:
	typedef T& reference;
	typedef const T& const_reference;
	explicit lockable(T &data, const boost::array<size_t,2> &shape)
	    : data_(data),
	      mutex_(new boost::mutex[shape[0]*shape[1]]),
	      array_(mutex_, shape) {}
	~lockable() { delete[] mutex_; }
	boost::mutex& lock(int i, int j) const {
	    return array_[j][i];
	}
	reference data(){return data_; }
	const_reference data()const {return data_; }
    private:
	T &data_;
	mutable boost::mutex *mutex_;
	mutable boost::multi_array_ref<boost::mutex,2> array_;
    };

    template<typename T>
    lockable<T> make_lockable(T &data) {
	return lockable <T>(data);
    }

    class fock_thread : boost::noncopyable {
    public:

	template<class B, class C, class S, class D, class F, class Q>
	static
	boost::thread* create(const B &blocks, const C &centers,
			      const S &screen,
			      const D &density, lockable<F> &fock,
			      const Parameters &parameters,
			      Q &queue) {
	    using boost::ref;
	    return new boost::thread
		(fock_thread::run<B, C, S, D, F, Q>,
		 blocks, ref(centers), ref(screen),
		 ref(density), ref(fock), parameters, ref(queue));
	}

	boost::thread thread_;
	template<class B, class C, class S, class D, class F, class Q>
	fock_thread(const B &blocks, const C &centers, const S &screen,
		    const D &density, lockable<F> &fock,
		    const Parameters &parameters, Q &queue)
	    : thread_(create(blocks, centers, screen, density, fock,
			     parameters, queue)) {}
	void join() { thread_.join(); }

    private:
	template<class B, class C, class S, class D, class F,
		 class Q = std::vector<boost::array<int,4> > >
	struct thread {
	private:
	    const B &blocks;
	    const C &centers;
	    const S &Screen;
	    const D &density;
	    lockable<F> &fock;
	    const Parameters &parameters;
	    Q quartets;
	    struct {
		typename F::matrix matrix;
		int index[2];
	    } fock_tls[6];
	public:
	    thread(const B &blocks, const C &centers, const S &Screen,
		   const D &density, lockable<F> &fock,
		   const Parameters &parameters)
		: blocks(blocks), centers(centers), Screen(Screen),
		  density(density), fock(fock), parameters(parameters) {}

	    template<class T>
	    void run(const T &task) {
	
	    BOOST_AUTO(ib, &blocks[task[0]]);
	    BOOST_AUTO(kb, &blocks[task[1]]);
	    BOOST_AUTO(jb, &blocks[task[2]]);
	    BOOST_AUTO(lb, &blocks[task[3]]);

	    const BOOST_TYPEOF(ib) block[] = { ib, jb, kb, lb };

	    BOOST_AUTO(screen, Screen(block, 1<<16));
	    quartets.clear();
	    screen(quartets);
	    if (quartets.empty()) return;

	    typedef adapter::rysq::Shell Shell;
	    typedef rysq::Quartet<rysq::Shell> Quartet;
	    Shell a(block[0]->shell());
	    Shell b(block[1]->shell());
	    Shell c(block[2]->shell());
	    Shell d(block[3]->shell());
	    Quartet quartet(a,b,c,d);

	    // if (quartet.L()/2 + 1 < 3) return;
	    
	    int base[] = { ib->start(), jb->start(), kb->start(), lb->start() };

	    rysq::hf::matrix_set<adapter::rysq::density_matrix> D_(base);
	    rysq::hf::matrix_set<adapter::rysq::fock_matrix> F_(base);

	    int k = 0;
	    foreach (const size_t (&index)[2],  rysq::index_list()) {
	        size_t i = index[0], j = index[1];
	        size_t im = std::distance(&blocks[0], block[i]);
	        size_t jm = std::distance(&blocks[0], block[j]);

	        D_.get(i,j) = adapter::rysq::density_matrix(density.m(im,jm));
	        D_.get(j,i) = adapter::rysq::density_matrix(density.m(jm,im));

		{
		    typename F::matrix &f = fock.data().m(im,jm);
		    fock_tls[k].matrix.resize(f.size1(), f.size2(),
		    			      f.block1(), f.block2());
		    fock_tls[k].matrix.clear();
		    fock_tls[k].index[0] = im;
		    fock_tls[k].index[1] = jm;
		    F_.get(i,j) = adapter::rysq::fock_matrix(fock_tls[k].matrix);
		    ++k;
		}
	    }

	    rysq::Fock::Parameters p(parameters.cscale,
				     parameters.xscale,
				     screen.cutoff());
	    rysq::Fock f(quartet);
	    while (!quartets.empty()) {
		f(centers, quartets, D_, F_, p);
		quartets.clear();
		screen(quartets);
	    }

	    for (int k = 0; k < 6; ++k) {
		int i = fock_tls[k].index[0], j = fock_tls[k].index[1];
		boost::lock_guard<boost::mutex> lock(fock.lock(i,j));
		fock.data().m(i,j) += fock_tls[k].matrix;
	    }

	}
	};

    public:
	template<class B, class C, class S, class D, class F, class Q>
	static void run(const B blocks, const C &centers, const S &screen,
			const D &density, lockable<F> &fock,
			const Parameters parameters, Q &queue) {

	    typename Q::value_type task;
	    timer t;
	    thread<B, C, S, D, F> thread(blocks, centers, screen,
					 density, fock, parameters);
	    while (true) {
	    	try { task = queue.backlog().pop(); }
	    	catch(...) {
	    	    try { task = queue.pop(); }
	    	    catch(...) { break; }
	    	}
	    	thread.run(task);
	    }
	    while(true) {
	    	try { task = queue.backlog().pop(); }
	    	catch(...) { break; }
		thread.run(task);
	    }

	    // std::cout << "CPU thread time: " << t << std::endl;
	}
    };


    struct fock_thread_group {
	static size_t max_threads() {
	    return std::max<size_t>(boost::thread::hardware_concurrency(), 1);
	}
	boost::thread_group group_;
	template<class B, class C, class S, class D, class F, class Q>
	fock_thread_group(const B &blocks, const C &centers, const S &screen,
			  const D &density, lockable<F> &fock,
			  const Parameters &parameters,
			  Q &queue, size_t threads = max_threads()) {
	    for (size_t i = 0; i < threads; ++i) {
		group_.add_thread(fock_thread::create
				  (blocks, centers, screen,
				   density, fock, parameters, queue));
	    }
	}
	void join_all() { group_.join_all(); }
    };


}

#endif /* HF_THREAD_HPP */
