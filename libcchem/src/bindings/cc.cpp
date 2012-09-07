#include "runtime.hpp"
#include "cc/cc.hpp"
#include "exception.hpp"

#include "bindings/bindings.h"
#include "bindings/bindings.hpp"


int CChem_cc_energy(int wf, double *E, const char* method) {
    const Wavefunction *pwf = bindings::any().find<Wavefunction*>(wf);
    CCHEM_ASSERT(pwf);
    try {
	cchem::cc::Energy e = cchem::cc::energy(*pwf, Runtime::rt(), method);
	E[0] = e["mp2"];
	E[1] = e["ccsd"];
	E[2] = e["ccsd[t]"];
	E[3] = e["ccsd(t)"];
    }
    catch (cchem::exception) { return 0; }
    return 1;
}
