#ifndef MP2_ENERGY_HPP
#define MP2_ENERGY_HPP

#include "core/wavefunction.hpp"
#include "core/integral.hpp"

namespace mp2 {
    double energy(const Wavefunction &wf, const Integral::Screening &screening);
}

#endif /* _AP2_MP2_HPP_ */
