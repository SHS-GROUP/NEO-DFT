#ifndef ARRAY_FORWARD_HPP
#define ARRAY_FORWARD_HPP

#include <boost/config.hpp>

namespace array {

    struct device_tag {};

    template<typename T, size_t N, class Tag = void, typename P = const T*>
    struct const_adapter;

    template<typename T, size_t N, class Tag = void>
    struct adapter;
}

#endif // ARRAY_FORWARD_HPP
