#include "multi_array/detail.hpp"
#include "multi_array/gpu/algorithm.hpp"

namespace multi_array {
namespace gpu {
namespace algorithm {
	    
    template void copy_traits<0123, const double, double, 4>::type::operator()();
    template void copy_traits<0231, const double, double, 4>::type::operator()();
	    
    template void pack::kernel<double, double>::type::operator()();

}	
}
}

