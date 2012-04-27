#!/usr/bin/python

import sys, re

#f = open("level1/copy.hpp").read()#sys.argv[1]).read()
f = open(sys.argv[1]).read()

f = re.sub("BLAS_LEVEL", "CUBLAS_LEVEL", f)
f = re.sub("namespace blas", "namespace cublas", f)

fail = "(throw std::runtime_error(std::string(__FILE__) + \":not implemented\"), 0);"
f = re.sub("// NOT FOUND\(\)", "%s;" % fail, f)
f = re.sub("(cublas(C|Z|Sc|Ic|Iz)[^;]*);", "%s;" % fail, f, re.M)
	   
#f = re.sub("cublasSc", "%s; //cublasSc" % fail, f)
#f = re.sub("#if defined BOOST_NUMERIC_BINDINGS_BLAS_CBLAS.*", "", f, re.M)

skip = False
blas = False
for line in f.split("\n"):
    if re.match("^namespace detail", line):
	line += "\nusing blas::detail::blas_option;"
	line += "\nusing blas::detail::default_order;"
    cblas_ = re.match(".*defined BOOST_NUMERIC_BINDINGS_BLAS_CBLAS", line)
    cublas_ = re.match(".*defined BOOST_NUMERIC_BINDINGS_BLAS_CUBLAS", line)
    else_ = re.match("#else", line)
    elif_ = re.match("#elif", line)
    endif_ = re.match("#endif", line)
    blas = any((blas, cblas_, cublas_)) and not endif_
    skip = any((skip, cublas_, cblas_, else_, elif_, endif_ and blas))

    if not skip:
	print re.sub("blas_option<", "blas::detail::blas_option<", line)
    else: print "////", line

    if re.match("#include <boost/numeric/bindings/blas/detail/cublas.h>", line):
	print "#include <boost/numeric/bindings/blas/detail/default_order.hpp>"
	print "#include <boost/numeric/bindings/cublas/exception.hpp>"
    if re.match(".*::invoke\(.*\);", line):
	print "    BOOST_NUMERIC_BINDINGS_CUBLAS_CHECK_STATUS();"

    
    skip = skip and not any((cublas_, elif_, endif_))


#print f

# for line in f:
#     line = re.sub("BLAS_LEVEL", "CUBLAS_LEVEL", line)
#     print line
