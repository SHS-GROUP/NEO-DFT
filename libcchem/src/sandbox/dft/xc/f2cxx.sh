#!/bin/sh

cat $1 | \
    sed \
    -e 's/(&/(/g' \
    -e 's/\, &/, /g' \
    -e 's/ \* /*/g' \
    -e 's/\*\*/**/g' \
    -e 's/ = \*/ = /g' \
    -e 's/(\*/(/g' \
    -e 's/\t\*/\t/g' \
    -e 's/    \*/    /g' \
    -e 's/doublereal */s/doublereal /g' \
    -e 's/pow_dd/pow/g' | \
    sed \
    -e 's/doublereal /const double /g' \
    -e 's/logical */const bool /g'

