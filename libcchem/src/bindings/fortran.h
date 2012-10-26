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

#define FORTRAN_FUNCTION(F, A, CF, CA) F ## _ A { return CF CA; }
#define FORTRAN_SUBROUTINE(F, A, CF, CA) void F ## _ A { CF CA; }

#define FORTRAN_CALL(F) F ## _

extern "C" {

#ifdef LIBCCHEM_WITH_INTEGER8
typedef int64_t integer_t; /**< @brief Use 64-bit integer */
#else 
typedef int32_t integer_t; /**< @brief Use 32-bit integer */
#endif



}
