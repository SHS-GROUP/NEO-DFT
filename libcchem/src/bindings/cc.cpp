#include "bindings/bindings.h"
#include "bindings/bindings.hpp"

#include "runtime.hpp"
#include "array/array.hpp"
#include "cc/cc.hpp"

void CC_sd_vvvv(int wf) {
    const Wavefunction *pwf = bindings::any().find<Wavefunction*>(wf);
    assert(pwf);

    const runtime &rt = runtime::rt();
    cc::Map<const cc::Array*> t;
	 
    t("i,a") = rt.get<Array<double>*>("cc.t(i,a)");
    t("a,b,i,j") = rt.get<Array<double>*>("cc.t(a,b,i,j)");
    t("i,j,a,b") = rt.get<Array<double>*>("cc.t(i,j,a,b)");
	 
    cc::sd::vvvv(*pwf, t,
		 *rt.get<Array<double>*>("cc.u(i,j,a,b)"),
		 *rt.get<Array<double>*>("cc.(n,n,o,o)"));
}

void CC_triples(size_t no, size_t nv,
		const double *eh, const double *ep,
		double* C) {
    const runtime &rt = runtime::rt();

    cc::Map<const Array<double>*> t, V;
    cc::Map<cc::Vector> e;

    t("i,a") = rt.get<Array<double>*>("cc.t(i,a)");
    // t("a,b,i,j") = rt.get<Array<double>*>("cc.t(a,b,i,j)");
    t("i,j,a,b") = rt.get<Array<double>*>("cc.t(i,j,a,b)");

    V("e,k,b,c") = rt.get<Array<double>*>("cc.V(e,k,b,c)");
    //V("a,b,c,i") = rt.get<Array<double>*>("cc.V(a,b,c,i)");
    V("j,k,i,a") = rt.get<Array<double>*>("cc.V(j,k,i,a)");
    V("i,j,a,b") = rt.get<Array<double>*>("cc.V(i,j,a,b)");

    e("i").resize(no);
    e("a").resize(nv);
    std::copy(eh, eh+no, e("i").begin());
    std::copy(ep, ep+nv, e("a").begin());

    cc::Triples::Correction c = cc::Triples()(no, nv, t, V, e);

    C[0] = c.ets;
    C[1] = c.etd;
    C[2] = c.ots;
    C[3] = c.otd;
}
