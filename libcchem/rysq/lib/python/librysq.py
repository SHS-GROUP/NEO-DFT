#!/usr/bin/env python

##
#####//file
#@brief Used for testing eris function
#

import math
import ctypes
import time
from ctypes import *

__librysq = ctypes.CDLL("librysq.so")

class Rysq_shell_t(Structure):
    _fields_ = [("type", c_int),
		("L", c_int),
		("min", c_int),
		("max", c_int),
		("size", c_int),
		("K", c_int),
		("nc", c_int),
		("a", c_void_p),
		("c", c_void_p)]

def Double(m, n = 1):
    return (ctypes.c_double*(m*n))()

def Rysq_type(a):
    if a == "s": return 1
    elif a == "p": return 2
    elif a == "sp": return 3
    elif a == "d": return 4
    else: return 2**(ord(a) - ord("f") + 3)


def Rysq_size(shell):
    if(isinstance(shell,str)):
        if(shell == "sp"):
            return Rysq_size("s") + Rysq_size("p")
        else:
            return Rysq_size(int(math.log(Rysq_type(shell),2)))
    else:
        L = shell
        return ((L + 1)*(L + 2))/2

def Rysq_ang_mom(l2):
    return int(math.log(l2,2))

def Rysq_shell(type, k = 1, a = None, C = None):
    global __librysq
    __librysq.Rysq_shell.restype = Rysq_shell_t
    if( type == 3 ):
	nc = 2
    else:
	nc = 1

    if not a: a = Double(k)
    
    if not C:
	C = ((ctypes.c_void_p)*2)()
	C[0] = cast(Double(k), c_void_p)
	C[1] = cast(Double(k), c_void_p)
	
    return  __librysq.Rysq_shell(type, k, a, nc, C);

def Rysq_shell_free(sh):
    global __librysq
    __librysq.Rysq_shell_free.restype = None
    __librysq.Rysq_shell_free.argtypes = [Rysq_shell_t]
    __librysq.Rysq_shell_free(sh)

        
def cuRysq_init(dev = 0):
    global __librysq
    __librysq.cuRysq_init.argtypes = [c_int]
    __librysq.cuRysq_init.restype = None
    __librysq.cuRysq_init(dev)

def cuRysq_eri(abcd,
	       flags = 0, tole = 0.0, k = 1,
	       ri=None, rj=None, rk=None, rl=None,
	       n = 1, index4 = None,
	       scale = 1.0, I=None):
              
    global __librysq
    __librysq.cuRysq_eri.argtypes = [c_int, c_double,
				     Rysq_shell_t, c_void_p,
				     Rysq_shell_t, c_void_p,
				     Rysq_shell_t, c_void_p,
				     Rysq_shell_t, c_void_p,
				     c_int, c_void_p, 
				     c_double, c_void_p]
    __librysq.cuRysq_eri.restype = c_int

    (bra, ket) = abcd.split("|")
    (a, b) = parseKet(bra)
    (c, d) = parseKet(ket)    
    
    type_a = Rysq_type(a)
    type_b = Rysq_type(b)
    type_c = Rysq_type(c)
    type_d = Rysq_type(d)

    if not ri : ri = Double(n*3)
    if not rj : rj = Double(3)
    if not rk : rk = Double(3)
    if not rl : rl = Double(3)

    erisize = Rysq_size(a)*Rysq_size(b)*Rysq_size(c)*Rysq_size(d)

    if not I : I = Double(n*erisize)

    sh_a = Rysq_shell(type_a, k = k)
    sh_b = Rysq_shell(type_b, k = k)
    sh_c = Rysq_shell(type_c, k = k)
    sh_d = Rysq_shell(type_d, k = k)

    if not index4: index4 = (ctypes.c_int*(n*4))()
    
    for ijkl in range(0, n):
	index4[0+ijkl] = ijkl
	index4[1+ijkl] = 0
	index4[2+ijkl] = 0
	index4[3+ijkl] = 0
		    
    sec_start = time.time()
    __librysq.cuRysq_eri(flags, tole,
			 sh_a, ri, sh_b, rj, sh_c, rk, sh_d, rl,
			 n, index4,
			 scale, I)
    sec_stop = time.time()
    sec_total = sec_stop - sec_start

    Rysq_shell_free(sh_a)
    Rysq_shell_free(sh_b)
    Rysq_shell_free(sh_c)
    Rysq_shell_free(sh_d)
    
    N = int((Rysq_ang_mom(type_a)+Rysq_ang_mom(type_b)+Rysq_ang_mom(type_c)+Rysq_ang_mom(type_d))/2) + 1
    flops = n*erisize*3*N
    GFLOPS = (flops/sec_total)/1000000000.0
    print "flops = %i" %(flops)
    print "time = %f sec" %(sec_total)
    print "GFLOPS = %f" %(GFLOPS)

def cuRysq_jk(abcd, ni, nj, nk, nl, 
	      flags = 0, tole = 0.0, scale = 1.0,
	      ri = None, rj = None, rk = None, rl = None,
	      Dij = None, Dkl = None, Dil = None, Djk = None,
	      Fij = None, Fkl = None, Fil = None, Fjk = None):

              
    global __librysq
    
    (bra, ket) = abcd.split("|")
    (a, b) = parseKet(bra)
    (c, d) = parseKet(ket)    
    
    type_a = Rysq_type(a)
    type_b = Rysq_type(b)
    type_c = Rysq_type(c)
    type_d = Rysq_type(d)

    if not ri : ri = Double(ni*3)
    if not rj : rj = Double(nj*3)
    if not rk : rk = Double(nk*3)
    if not rl : rl = Double(nl*3)

    erisize = Rysq_size(a)*Rysq_size(b)*Rysq_size(c)*Rysq_size(d)
    nij = Rysq_size(a)*Rysq_size(b)
    nkl = Rysq_size(c)*Rysq_size(d)
    nil = Rysq_size(a)*Rysq_size(d)
    njk = Rysq_size(b)*Rysq_size(c)

    if not Dij : Dij = Double(nij*ni*nj)
    if not Dkl : Dkl = Double(nkl*nk*nl)
    if not Dil : Dil = Double(nil*ni*nl)
    if not Djk : Djk = Double(njk*nj*nk)
    if not Fij : Fij = Double(nij*ni*nj)
    if not Fkl : Fkl = Double(nkl*nk*nl)
    if not Fil : Fil = Double(nil*ni*nl)
    if not Fjk : Fjk = Double(njk*nj*nk)

    sec_start = time.time()
    __librysq.cuRysq_jk(flags, c_double(tole),
			Rysq_shell(type_a), ni, ri,
			Rysq_shell(type_b), nj, rj,
			Rysq_shell(type_c), nk, rk,
			Rysq_shell(type_d), nl, rl,
			Dij, Dkl, Dil, Djk,
			c_double(1.0), Fij, Fkl, Fil, Fjk)
    sec_stop = time.time()
    sec_total = sec_stop - sec_start

    N = int((Rysq_ang_mom(type_a)+Rysq_ang_mom(type_b)+Rysq_ang_mom(type_c)+Rysq_ang_mom(type_d))/2) + 1
    flops = ni*nj*nk*nl*(erisize*3*N+nij*nkl*2*4)
    GFLOPS = (flops/sec_total)/1000000000.0
    print "flops = %i" %(flops)
    print "time = %f sec" %(sec_total)
    print "GFLOPS = %f" %(GFLOPS)

def Rysq_init():
    global __librysq
    __librysq.Rysq_init.argtypes = None
    __librysq.Rysq_init.restype = c_int
    __librysq.Rysq_init()

def Rysq_finalize():
    global __librysq
    __librysq.Rysq_finalize.argtypes = None
    __librysq.Rysq_finalize.restype = c_int
    __librysq.Rysq_finalize()

def Rysq_eri(abcd,
	     flags = 0, tole = 0.0, k = 1,
	     ri = None, rj = None, rk = None, rl = None,
	     n = 1, index4 = None,
	     scale = 1.0, I = None):
              
    global __librysq
    __librysq.Rysq_eri.argtypes = [c_int, c_double,
			     Rysq_shell_t, c_void_p,
			     Rysq_shell_t, c_void_p,
			     Rysq_shell_t, c_void_p,
			     Rysq_shell_t, c_void_p,
			     c_int, c_void_p, 
			     c_double, c_void_p]
    __librysq.Rysq_eri.restype = c_int

    (bra, ket) = abcd.split("|")
    (a, b) = parseKet(bra)
    (c, d) = parseKet(ket)    
    
    type_a = Rysq_type(a)
    type_b = Rysq_type(b)
    type_c = Rysq_type(c)
    type_d = Rysq_type(d)

    if not ri : ri = Double(n*3)
    if not rj : rj = Double(3)
    if not rk : rk = Double(3)
    if not rl : rl = Double(3)

    erisize = Rysq_size(a)*Rysq_size(b)*Rysq_size(c)*Rysq_size(d)

    if not I : I = Double(n*erisize)

    sh_a = Rysq_shell(type_a, k = k)
    sh_b = Rysq_shell(type_b, k = k)
    sh_c = Rysq_shell(type_c, k = k)
    sh_d = Rysq_shell(type_d, k = k)

    if not index4: index4 = (ctypes.c_int*(n*4))()
    
    for ijkl in range(0, n):
	index4[0+ijkl] = ijkl
	index4[1+ijkl] = 0
	index4[2+ijkl] = 0
	index4[3+ijkl] = 0
		    
    sec_start = time.time()
    __librysq.Rysq_eri(flags, tole,
 		     sh_a, ri, sh_b, rj, sh_c, rk, sh_d, rl,
 		     n, index4,
 		     scale, I)
    sec_stop = time.time()
    sec_total = sec_stop - sec_start

    Rysq_shell_free(sh_a)
    Rysq_shell_free(sh_b)
    Rysq_shell_free(sh_c)
    Rysq_shell_free(sh_d)


    N = int((Rysq_ang_mom(type_a)+Rysq_ang_mom(type_b)+Rysq_ang_mom(type_c)+Rysq_ang_mom(type_d))/2) + 1
    flops = n*erisize*3*N*(k*k*k*k)
    GFLOPS = (flops/sec_total)/1000000000.0
    print "flops = %i" %(flops)
    print "time = %f sec" %(sec_total)
    print "GFLOPS = %f" %(GFLOPS)


def parseKet(ket):
    if len(ket) == 2: return (ket[0], ket[1])
    elif len(ket) == 3:
	if ket.startswith("sp"): return ("sp", ket[2])
	else: return (ket[0], "sp")
    else: return ("sp", "sp")
    
    
           
import sys
def main():
    abcd = sys.argv[1]
    num_centers = int(sys.argv[2])
    if len(sys.argv) == 4: k = int(sys.argv[3])
    if len(sys.argv) == 7: device = int(sys.argv[6])
    else: device = 0
    cuRysq_init(dev = device)
    cuRysq_eri(abcd, n = num_centers)
#     Rysq_init()
#     Rysq_eri(abcd, n = num_centers, k = k)
#     Rysq_finalize()

if __name__ == "__main__":
    main()
            
