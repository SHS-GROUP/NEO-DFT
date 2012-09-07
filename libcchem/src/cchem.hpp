#ifndef CCHEM_HPP
#define CCHEM_HPP

#include "runtime.hpp"
#include "integrals/integrals.hpp"

namespace cchem {

    inline void initialize() {
	Runtime::rt();
	integrals::initialize();
    }

}

#endif /* CCHEM_HPP */

