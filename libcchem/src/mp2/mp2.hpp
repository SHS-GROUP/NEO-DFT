#ifndef MP2_MP2_HPP
#define MP2_MP2_HPP

#include "runtime.hpp"
#include "core/wavefunction.hpp"

namespace cchem {
namespace mp2 {

    double energy(Wavefunction wf, Runtime &rt);

}
} // namespace cchem

#endif /* MP2_MP2_HPP */
