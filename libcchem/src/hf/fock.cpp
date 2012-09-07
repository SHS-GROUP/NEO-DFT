/**
 * @file   
 * @brief  Fock matrix operator implementation
 */


#include "hf/work_queue.hpp"
#include "hf/hf.hpp"
#include "hf/thread.hpp"
#include "hf/device_thread.hpp"

#include "runtime.hpp"
#include "omp.hpp"

#include "matrix/meta.hpp"

#include <algorithm>
#include <boost/typeof/typeof.hpp>
#include <boost/utility.hpp>
#include <boost/math/special_functions/pow.hpp>
#include "boost/utility/profiler.hpp"


template
void hf::fock<hf::Matrix>(const Basis &basis, 
			  const matrix::meta_matrix<BlockMatrix> &D,
			  matrix::meta_matrix<BlockMatrix> &F,
			  const Parameters &parameters,
			  Parallel &pe,
			  const Matrix *Kmax, const Matrix *Dmax);

template<class M>
void hf::fock(const Basis &basis, 
	      const matrix::meta_matrix<BlockMatrix> &D,
	      matrix::meta_matrix<BlockMatrix> &F,
	      const Parameters &parameters,
	      Parallel &pe,
	      const M *Kmax, const M *Dmax) {

    timer t;

    Screen<M> screen(*Kmax, *Dmax, parameters.cutoff);

    rysq::initialize();

    pe.task().reset();
    work_queue_lock queue_lock;
    work_queue queue(queue_lock, basis.blocks().size(), pe.task());
    lockable<BOOST_TYPEOF(F)> lockable_fock(F, F.meta().shape());
    std::vector<Basis::Center> centers(basis.centers());

    Runtime &rt = Runtime::rt();

#pragma omp parallel
    {

	hf::Thread<BlockMatrix> thread;
	hf::Thread<BlockMatrix>::Device *device = 0;

#pragma omp master
	if (pe.rank() == 0) {
	    // rt.cout() << "Threads: " << omp::num_threads() << std::endl;
	    BOOST_PROFILE_REGISTER_THREAD;
	}

#if HF_CUDA
	{
	    std::vector<int> devices = rt.devices(pe.node());
	    if (omp::thread() < (int)devices.size()) {
		// ::cuda::thread_exit();
		try {
		    cuda::set_device(devices.at(omp::thread()));
		    device = new hf::Thread<BlockMatrix>::Device();
		}
		catch (::cuda::error &e) {
		    CCHEM_MESSAGE("%s\n", e.what());
		    delete device;
		    device = NULL;
		}
	    }
	    if (device) {
		BOOST_PROFILE_LINE;
		(*device)(basis.blocks(), centers, D, lockable_fock,
			  screen, queue, parameters);
	    }
	}
#endif

	if (!device) {
	    BOOST_PROFILE_LINE;
	    thread(basis.blocks(), centers, D, lockable_fock,
		   screen, queue, parameters);
	}

#if HF_CUDA
	delete device;
#endif
	
    }

    {
	BOOST_AUTO(cout, rt.cout());
	BOOST_PROFILE_DUMP(cout);
    }

    pe.barrier();

    // sum fock matrix across processors
    for (size_t j = 0; j < F.array().size(); ++j) {
    	for (size_t i = 0; i < F.array()[0].size(); ++i) {
    	    BOOST_AUTO(&f, F.m(i,j));
    	    pe.reduce("+", f.data(), f.size1()*f.size2());
    	}
    }

    matrix::symmeterize(F); // symmetrize fock matrix

    return;
}
