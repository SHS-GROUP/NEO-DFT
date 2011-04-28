#ifndef DDI_UTIL_H
#define DDI_UTIL_H

/**
   @file
   @brief Miscellaneous utility functions.
*/

#include <stdlib.h>

/**
   @brief Join two strings by a delimiter.
   @note The strings are joined by the delimiter only if the length of the first string
   is greater than 0.
   @param[in,out] to First and destination string.
   @param[in] from Second string.
   @param[in] delim Delimiter string.
*/
void join(char *to, const char *from, const char *delim);

/** 
    @brief Converts at most len-1 characters of Fortran string to NULL terminated C string.
    @param[out] cstr C string.
    @param[in] fstr Fortran string.
    @param[in] len String length including NULL character.
*/
void strf2c(char *cstr, char *fstr, size_t len);

#endif /* DDI_UTIL_H */
