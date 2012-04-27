/**
@file
@brief librysq fortran bindings
@details The fortran binding provides lower case subroutines with and without trailing underscore for compatibility with different fortran compilers.
The default size of fortran integers is 4 bytes. If LIBCCHEM_WITH_INTEGER8 is defined, the size is 8 bytes. 
*/

#ifdef HAVE_CONFIG_H 
#include "config.h"
#endif

#include <stdint.h>

extern "C" {

#ifdef LIBCCHEM_WITH_INTEGER8
typedef int64_t Integer; /**< @brief Use 64-bit integer */
#else 
typedef int32_t Integer; /**< @brief Use 32-bit integer */
#endif

}
