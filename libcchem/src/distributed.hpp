#ifndef _DISTRIBUTED_HPP_
#define _DISTRIBUTED_HPP_

#include "distributed/forward.hpp"
#include "distributed/array.hpp"
#include "distributed/matrix.hpp"

template<typename T = int>
struct index_range : std:: vector<T> {
    typedef index_range self;
    struct range {
	T first, last, inc;
	range(T idx) : first(idx), last(first), inc(1) {}
    };
    self& operator()(T idx) { this->push_pack(range(idx)); return *this; }
};

#endif /* _DISTRIBUTED_HPP_ */
