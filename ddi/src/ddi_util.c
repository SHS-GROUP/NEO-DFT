#include "ddi_util.h"
#include <string.h>
#include <stdlib.h>

/** @see ddi_util.h */
void join(char *to, const char *from, const char *delim) {
    if (strlen(to) > 0) strcat(to,delim);
    strcat(to,from);
}

/** @see ddi_util.h */
void strf2c(char *cstr, char *fstr, size_t len) {
    strncpy(cstr, fstr, len-1);
    cstr[len-1] = '\0';
}
