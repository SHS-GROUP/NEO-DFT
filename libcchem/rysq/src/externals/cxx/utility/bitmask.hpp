#ifndef _CXX_UTILITY_BITMASK_HPP_
#define _CXX_UTILITY_BITMASK_HPP_

#include <stdlib.h>
#include <stdint.h>

namespace cxx {
    namespace utility {

	template<typename T>
	int bitmask(const T &b0, const T &b1, const T &b2, const T &b3) {
	    return (bool(b0) | bool(b1) << 1 | bool(b2) << 2 | bool(b3) << 3);
	}

	inline int32_t set(uint8_t b, size_t p) { return (int32_t(b) << 8*p); }

	inline int32_t pack(int8_t b0, int8_t b1, int8_t b2, int8_t b3) {
	    return set(b0,0) | set(b1,1) | set(b2,2) | set(b3,3);
	}

    }
}

#endif /* _CXX_UTILITY_BITMASK_HPP_ */
