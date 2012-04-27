#ifdef HAVE_CONFIG_H 
#include "config.h"
#endif

#include <stdint.h>

extern "C" {

#ifdef RYSQ_WITH_INTEGER8
    typedef int64_t Integer; /**< @brief Use 64-bit integer */
#else 
    typedef int32_t Integer; /**< @brief Use 32-bit integer */
#endif

}
