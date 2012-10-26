#ifndef BOOST_CUDA_FORWARD_HPP
#define BOOST_CUDA_FORWARD_HPP

namespace boost {
namespace cuda {

    struct thread;

    struct flags {
	enum flag {
	    none = 0,
	    schedule_auto = 1,
	    schedule_spin = 2,
	    schedule_yield = 4,
	    blocking_sync = 8,
	    map_host = 16
	};
	static void set(flag value);
    };


    struct cache {
	enum config {
	    none = 0,
	    l1 = 1,
	    shared = 2
	};
	static void set(config value);
    };

}
}

#endif /* BOOST_CUDA_FORWARD_HPP */
