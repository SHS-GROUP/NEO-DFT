#!/usr/bin/env sage

import sys
import operator
import pprint
from sympy import Symbol, expand, sympify, symbols, Integer
from sympy.simplify.cse_main import numbered_symbols
from shell import *
import sage.all

def mathematica_to_sympy(e):
    return sympify(str(e).replace("^", "**").replace("\n", ""))



def simplify(Q, subs):
    assum = ", ".join([ "%s -> %s" % (v,k) for (k,v) in subs ])
    assum = assum.replace("**", "^")
    math = sage.all.Mathematica(maxread = 1000)
    P = []
    for (q,e) in Q:
	e = str(e).replace("**", "^")
	#cmd = "InputForm[FullSimplify[FullSimplify[%s]] /. { %s }]" % (e, assum)
	cmd = "InputForm[FullSimplify[%s, {%s}] /. { %s }]" % (e, assum, assum)
	#print cmd, math.eval(cmd)
	P.append((q,mathematica_to_sympy(math.eval(cmd))))
    return P

def factor(Q):
    math = sage.all.Mathematica(maxread = 1000)
    P = []
    for (q,e) in Q:
	e = str(e).replace("**", "^")
	cmd = "InputForm[Factor[%s]]" % (e)
	#print math.eval(cmd)
	P.append((q,mathematica_to_sympy(math.eval(cmd))))
    return P

def sympy_to_cxx(t):
    if not t: return ""
    if t.is_Atom: return str(t)
    ts = map(sympy_to_cxx, t.args)
    if t.is_Mul: return "*".join(ts)
    if t.is_Add: return "(%s)" % " + ".join(ts)
    if t.is_Pow: return "pow(%s,%s)" % (ts[0], ts[1])
    raise "sympy expression not implemented"

def tokens(p):
    if p.is_Add: return p.args
    return [t for t in p.args if t.is_Add]

def optimize(exprs, subs, rearrange = False):
    # map expressions with tokens to token
    token_exprs = {}
    for (q,e) in exprs:
	for t in tokens(e):
	    if not t in token_exprs: token_exprs[t] = []
	    token_exprs[t] += [e]

    # expressions with common tokens/terms
    T = numbered_symbols("f")
    subs += [(T.next(), k) for (k,v) in token_exprs.items() if len(v) > 1]
    
    # replace common terms with Tn symbol
    subs_exprs = []
    for (q,e) in exprs:
	for (T,t) in subs: e = e._eval_subs(t, T)
	subs_exprs += [(q,e)]

    if rearrange:
	cmp = lambda x,y: len(x[1].atoms(Symbol)) - len(y[1].atoms(Symbol))
	Q = sorted(subs_exprs, cmp, reverse = True)
	subs_exprs = [Q.pop()]
	while Q:
	    a = subs_exprs[-1][1].atoms(Symbol)
	    (p,s) = Q[-1]
	    for (q,t) in Q[:-1]:
		if len(t.atoms(Symbol) | a) < len(s.atoms(Symbol) | a): (p,s) = (q,t)
	    subs_exprs.append((p,s))
	    Q.remove((p,s))
	subs_exprs.reverse()

    return (subs_exprs, subs)

def resolve_symbols(E, S, T):
    R = []
    for e in E.atoms(Symbol):
	if e in T: R += [(e, None)]
	else:
	    s = dict(S)[e]
	    R += [(e, s)] + resolve_symbols(s, S, T)
    return R
    
class Quadrature:
    (B00, B10, B01) = symbols("B00", "B10", "B01")
    (Cx, Dx) = symbols("Cx", "Dx")
    (Cy, Dy) = symbols("Cy", "Dy")
    (Cz, Dz) = symbols("Cz", "Dz")
    (xij, yij, zij) = symbols("xij", "yij", "zij")
    (xkl, ykl, zkl) = symbols("xkl", "ykl", "zkl")

    primary = { B00: None, B10: None, B01: None,
		Cx: None, Cy: None, Cz: None,
		Dx: None, Dy: None, Dz: None,
		xij: None, yij: None, zij: None,
		xkl: None, ykl: None, zkl: None }
    
    (Px, Qx, Rx, Sx, Tx) = symbols("Px", "Qx", "Rx", "Sx", "Tx")
    (Py, Qy, Ry, Sy, Ty) = symbols("Py", "Qy", "Ry", "Sy", "Ty")
    (Pz, Qz, Rz, Sz, Tz) = symbols("Pz", "Qz", "Rz", "Sz", "Tz")

    (Ix, Kx) = symbols("Ix", "Kx")
    (Iy, Ky) = symbols("Iy", "Ky")
    (Iz, Kz) = symbols("Iz", "Kz")

    subs = { Px : B10 + Cx**2, Qx : B00 + Cx*Dx, Rx : B01 + Dx**2,
	     Py : B10 + Cy**2, Qy : B00 + Cy*Dy, Ry : B01 + Dy**2,
	     Pz : B10 + Cz**2, Qz : B00 + Cz*Dz, Rz : B01 + Dz**2 }

    hints = { Ix : Cx + Symbol("xij"), Kx : Dx + Symbol("xkl"),
	      Iy : Cy + Symbol("yij"), Ky : Dy + Symbol("ykl"),
	      Iz : Cz + Symbol("zij"), Kz : Dz + Symbol("zkl"),
	      Sx : B10 + Cx*Ix, Tx : B01 + Dx*Kx,
	      Sy : B10 + Cy*Iy, Ty : B01 + Dy*Ky,
	      Sz : B10 + Cz*Iz, Tz : B01 + Dz*Kz }

    
    @staticmethod
    def recurrence(i, j, q):
	(B00, B10, B01) = symbols("B00", "B10", "B01")
	C = Symbol("C" + q)
	D = Symbol("D" + q)
	P = Symbol("P" + q)
	Q = Symbol("Q" + q)
	R = Symbol("R" + q)
	G = Quadrature.recurrence 
	if i < 0 or j < 0: return 0
	if (i,j) == (0,0): return Integer(1)
# 	if (i,j) == (2,0): return B10 + C**2 #P
# 	if (i,j) == (1,1): return Q
# 	if (i,j) == (0,2): return R
	if j == 0: return (i-1)*B10*G(i-2,0,q) + C*G(i-1,0,q)
	return (j-1)*B01*G(i,j-2,q) + D*G(i,j-1,q) + i*B00*G(i-1,j-1,q)

    @staticmethod
    def transfer(i, j, k, l, q):
	qij = Symbol(q + "ij")
	qkl = Symbol(q + "kl")
	transfer = Quadrature.transfer
	if j == 0 and l == 0: return Quadrature.recurrence(i,k,q)
	if j > 0: return qij*transfer(i,j-1,k,l,q) + transfer(i+1,j-1,k,l,q)
	if l > 0: return qkl*transfer(i,j,k,l-1,q) + transfer(i,j,k+1,l-1,q)

    @staticmethod
    def expr1(q):
	transfer = Quadrature.transfer
	xs = [o[0] for o in q] + ["x"]
	ys = [o[1] for o in q] + ["y"]
	zs = [o[2] for o in q] + ["z"]
	return transfer(*xs)*transfer(*ys)*transfer(*zs)

    @staticmethod
    def exprs(quartets):
	return [(q, Quadrature.expr1(q)) for q in quartets]

    @staticmethod
    def generate1(a, b = "s", c = "s", d = "s"):
	if type(a) is str: a = Shell(a)
	if type(b) is str: b = Shell(b)
	if type(c) is str: c = Shell(c)
	if type(d) is str: d = Shell(d)

	Q = Quadrature.exprs(tuple(product(a,b,c,d)))
	Q = factor(Q)
	(Q,S) = optimize(Q, [], rearrange = True)

	# simplify common factors
	subs =  Quadrature.hints.items() + Quadrature.subs.items()
	S = simplify(S, subs)
	    
	# remove identity elements
	for (T,t) in S:
	    if len(t.args) == 0: Q = [(q, e._eval_subs(T, t))for (q,e) in Q]
	S = [(T, t)for (T, t)in S if len(t.args) > 0]

	# simplify expressions
	S += Quadrature.subs.items() + Quadrature.hints.items()
	Q = simplify(Q, S)

	resolved = set()
	for (q,e) in Q:
	    resolved.update(resolve_symbols(e, S, Quadrature.primary))

	return (Q, tuple(resolved))


def cxx(Q, S):
    cxx_exprs = tuple([(k, (sympy_to_cxx(e), map(repr, e.atoms(Symbol))))
		      for (k,e) in Q])
    cxx_subs = dict([(repr(s), sympy_to_cxx(t)) for (s,t) in S])
    return (cxx_subs, cxx_exprs)

def main():
    S =  map(str, Shell.range(last = 3, sp = True))
    Q = {}
    for q in product(S,S,S,S):
	L = sum(map(Shell.ang_mom, q));
 	(a, b, c, d) = map(Shell.size, q)
	size = a*b*c*d
	#(a, b, c, d) = map(Shell.value, q)
	conditions = []
	conditions += [ L > 4 and size > 160 ]
	conditions += [ b > a, d > c, c > a ]
	conditions += [ a == c and b < d ]
	# for gpu, provided non-permuted small kernels
	if any(conditions) and size > 40: continue
	#print q
	sys.stderr.write(str(q) + "\n")
	Q[q] = cxx(*(Quadrature.generate1(*q)))

    print "# quadrature expressions: Q[ShellQuartet] = ( {Symbol:Expressions},"
    print "# ((i, j, k, l), (Expression, [Symbols])))"
    pprint.pprint(Q)

if __name__ == "__main__":
    main()

	


    

