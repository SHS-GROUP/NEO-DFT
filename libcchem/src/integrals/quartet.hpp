#ifndef INTEGRALS_QUARTET_HPP
#define INTEGRALS_QUARTET_HPP

#include "basis/basis.hpp"
#include "foreach.hpp"

namespace integrals {


    template<typename T>
    struct Quartet {
	T p, q, r, s;
    };

    template<>
    struct Quartet<const Basis::Shell::Data&> {
	const Basis::Shell::Data &p, &q, &r, &s;
	size_t size() const {
	    return (p.size()*q.size()*r.size()*s.size());
	}
    };


}

#endif // INTEGRALS_QUARTET_HPP
