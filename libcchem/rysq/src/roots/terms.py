#!/usr/bin/python2.6

import sys
import ast
import sympy

def tokenize(expression):
    if isinstance(expression, sympy.Add):
	tokens = expression.as_two_terms()
	return  tokenize(tokens[0]) + tokenize(tokens[1])
    else: return [ expression ]

def to_list(terms):
    l = []
    for (k,v) in terms.items():
	l += ["0"]*(max(0, k - len(l) + 1))
	l[k] = str(v).replace("*em", "E-")
    return l

def token(string):
    terms = sympy.sympify(string).expand()
    terms = tokenize(terms)

    inverse = {}
    polynomial = {}
    for term in terms:
	ce = term.as_coeff_exponent(sympy.Symbol("x"))
	c = ce[0]
	e = int(ce[1])
	if e < 0: inverse[-e] = c
	else:
	    if e in polynomial: polynomial[e] += c
	    else: polynomial[e] = c

    return (to_list(polynomial), to_list(inverse))



lines = "".join([line.strip("\n ") for line in sys.stdin.readlines()])
strings = lines.split(";")

polynomials = []
inverses = []

for string in strings:
    string = string.replace("\t", "").strip().lower()
    string = string.replace("y", "x")
    string = string.replace("e-", "*em")
    string = string.partition("=")[2]
    if string:
	(polynomial, inverse) = token(string)
	polynomials += [polynomial]
	inverses += [inverse]

strings = []
m = max(len(l)for l in polynomials)
for l in polynomials:
    l += ["0.0"]*(m - len(l))
    strings += [ "{ " + ", ".join(l[0:])+ " }" ]

declaration = "const double CP[" + str(len(strings)) + "][" + str(m) + "] = "
print declaration + "{" + ",\n".join(strings) + "};"

strings = []
m = max(len(l)for l in inverses)
for l in inverses:
    l += ["0.0"]*(m - len(l))
    strings += [ "{ " + ", ".join(l[1:])+ " }" ]

# declaration = "double CR[" + str(len(strings)) + "][" + str(m-1) + "] = "
# print declaration + "{" + ",\n".join(strings) + "};"

