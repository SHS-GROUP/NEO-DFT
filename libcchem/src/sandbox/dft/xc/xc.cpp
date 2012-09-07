#include "dft/xc/functional.hpp"
#include "dft/xc/library.hpp"

#include <string>
#include <stdexcept>

#include <boost/fusion/include/make_vector.hpp>

namespace dft {
namespace xc {

    functional* functional::new_(const std::string &name) {

#define FUNCTIONAL(NAME, F)				\
	if (name == NAME)				\
	    return new_(boost::fusion::make_vector F);

	FUNCTIONAL("b3lyp", (slater(0.08), b88(0.72), lyp(0.81), vwn5(0.19)));

	throw std::runtime_error(std::string(__FILE__) + ": " +
				 "invalid functional \"" + name +  "\"");
    }

}
}

