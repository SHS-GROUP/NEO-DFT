#ifndef HF_FOCK_HPP
#define HF_FOCK_HPP

/**
 * @file   
 * @brief  Fock matrix operator implementation
 */

#include <algorithm>

#include "hf/work_queue.hpp"
#include "hf/hf.hpp"
#include "hf/thread.hpp"
#include "runtime.hpp"

#include "util/numeric.hpp"
#include "matrix/meta.hpp"

#include "hf/cuda.hpp"

#include <boost/typeof/typeof.hpp>
#include <boost/utility.hpp>
#include <boost/math/special_functions/pow.hpp>
#include "boost/utility/profiler.hpp"

#ifdef HAVE_GA
#include "distributed.hpp"
#endif


template
void hf::fock<Matrix>(const Basis &basis, 
		      const matrix::meta_matrix<BlockMatrix> &D,
		      matrix::meta_matrix<BlockMatrix> &F,
		      const Parameters &parameters,
		      parallel::counter<size_t> &counter,
		      const Matrix *Kmax, const Matrix *Dmax);

template<class M>
void hf::fock(const Basis &basis, 
	      const matrix::meta_matrix<BlockMatrix> &D,
	      matrix::meta_matrix<BlockMatrix> &F,
	      const Parameters &parameters,
	      parallel::counter<size_t> &counter,
	      const M *Kmax, const M *Dmax) {

    timer t;

    Screen<M> screen(*Kmax, *Dmax, parameters.cutoff);

    rysq::initialize();

    work_queue_lock queue_lock;
    work_queue queue(queue_lock, basis.blocks().size(), counter);
    lockable<BOOST_TYPEOF(F)> lockable_fock(F, F.meta().shape());

    std::vector<Basis::Center> centers(basis.centers());

    runtime rt("HF");
    int num_threads = rt.get<int>("threads");
    std::vector<int> devices = rt.get< std::vector<int> >("gpu::devices");
    num_threads -= (devices.size() > 0);
    num_threads = std::max(num_threads,1);

    fock_thread_group host(basis.blocks().begin(), centers, screen,
			   D, lockable_fock, parameters,
			   queue, num_threads);

#if HF_CUDA
    if (!devices.empty()) {
	timer t;
    	cuda::fock_thread_group cuda(basis.blocks().begin(), centers, screen,
    				     D, lockable_fock, parameters,
				     queue, devices);
    	cuda.join_all();

	if (rt.profile()) {
	    std::ostream_iterator<int> it(std::cout, ",");
	    std::cout << "gpu time (devices: ";
	    std::copy(devices.begin(), devices.end()-1, it);
	    std::cout << devices.back();
	    std::cout << "): " << t << std::endl;
	}

    }
#endif
    host.join_all();

    if (!queue.backlog().empty()) {
	fock_thread_group host(basis.blocks().begin(), centers, screen,
			       D, lockable_fock, parameters,
			       queue, num_threads);
	host.join_all();
    }

    if (rt.profile()) {
	std::cout << "hf time (" << num_threads << " threads): "
		  << t << std::endl;
    }

    //std::cout << boost::utility::global_profiler() << std::endl;
    boost::utility::global_profiler().clear();

    matrix::symmeterize(F); // symmetrize fock matrix

    // surface.render(K);
    //matrix::plot::surface(F);

    return;
}


#endif /* HF_FOCK_HPP */

