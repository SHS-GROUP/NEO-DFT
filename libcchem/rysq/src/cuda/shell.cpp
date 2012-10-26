#include "rysq-cuda.hpp"

#include <boost/noncopyable.hpp>
#include "boost/utility/profiler.hpp"

#include <vector>
#include <map>
#include "cuda.hpp"

namespace rysq {
namespace cuda {

    Shell::Shell(const rysq::Shell &shell) {
	// BOOST_PROFILE_LINE;

	std::vector<Primitive> host;

	for (int i = 0; i < shell.K; ++i) {
	    Primitive p = { shell(i), shell(i,0) };
	    if (shell == rysq::SP) {
		p.C = shell(i,1);
		p.Cs = shell(i,0)/p.C;
	    }
	    host.push_back(p);
	}

	// // pack exponents
	// for (int i = 0; i < shell.K; ++i) {
	//     host.push_back(shell(i));
	// }
	// // pack coefficients
	// for (int j = 0; j < shell.nc; ++j) {
	//     for (int i = 0; i < shell.K; ++i) {
	// 	host.push_back(shell(i,j));
	//     }
	// }

	namespace rt = runtime;
	Primitive* data = static_cast<Primitive*>
	    (rt::malloc(host.size()*sizeof(Primitive)));
	rt::copy(host.size(), &host[0], data, rt::host_to_device);
	this->data_ = data;
	this->type_ = shell.type;
	this->K_ = shell.K;
	// std::cout << "shell: " << data_ << std::endl;
    }

    Shell::~Shell() { runtime::free(const_cast<Primitive*>(this->data_)); }

}
}
