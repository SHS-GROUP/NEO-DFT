#!/usr/bin/python

## @file
# Shell classes

from math import sqrt


def product(*args, **kwds):
    # almost identical to python product function, except for ordering
    # product('ABCD', 'xy') --> Ax Bx Cx Dx Ay By Cy Dy
    # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
    pools = map(tuple, args) * kwds.get('repeat', 1)
    result = [[]]
    for pool in pools:
        result = [x+[y] for y in pool for x in result]
    for prod in result:
        yield tuple(prod)

def coefficient_index2(i, A, B, C, D):
    m = len(A)*len(B)
    ik = int(A.hybrid and (i%m)%len(A) > 0)
    jk = int(B.hybrid and (i%m) >= len(A))
    kk = int(C.hybrid and (i/m)%len(C) > 0)
    lk = int(D.hybrid and (i/m) >= len(C))
    return [ ik + jk*(A.nc), kk + lk*(C.nc) ]

def coefficient_index(i, A, B, C, D):
    (ij, kl) = coefficient_index2(i, A, B, C, D)
    return ij + kl*(A.nc*B.nc)


## Shell class, encapsulates information about shells and orbitals
# The shell class is derived from tuple, the index operator [] allows to access individual orbitals of a shell.
class Shell(tuple):
	"""
	Shell container
	Shell attributes;
	type 
	"""
	
	## Returns true if the shell type is SP
	# @param type shell type
	def is_sp(type):
		"""
		Check to see if the Shell is sp
		"""
		if(str(type) == "sp"): return True
		return False

	## Returns shell type corresponding to angular momentum
	# @param L angular momentum
	def type(L):
		"""
		Returns the shell type
		"""
		if L == 0: return "s"
		elif L == 1: return "p"
		elif L == 2: return "d"
		else: return str(chr(ord('f') + L - 3))

	## Returns angular momentum corresponding to a shell type
	# @param a Shell type
	def ang_mom(a):
		"""
		Returns the angular momentum of a shell
		"""
		if a == "s": return 0
		elif a == "p": return 1
		elif a == "sp": return 1
		elif a == "d": return 2
		else: return ord(a) - ord("f") + 3

	@staticmethod
	def value(shell):
	    shell = Shell(shell)
	    return ((1 << shell.L) | (shell.hybrid)*((1 << shell.L) - 1)) - 1
    

	
	## 	Returns the shell size
	# @param type Shell type
	def size(type):
		"""
		Return the shell size
		"""
		L = Shell.ang_mom(type)
		if Shell.is_sp(type): return Shell.size("s")+Shell.size("p")
		return ((L+1)*(L+2))/2

	## Returns a range of shell objects between angular momentums from first to last inclusive
	# @param first angular momentum
	# @param last angular momentum
	# @param sp Include SP shell
	def range(first = 0, last = 0, sp = False):
		"""
		Returns all the shells up to a given L
		"""
		types = [Shell.type(L) for L in range(first, last+1)]
		if sp: types.append(Shell.type(0)+Shell.type(1))
		return tuple([Shell(t) for t in types])

	## Returns 2-tuple of shell types in the given bra
	# @param bra Bra pair <ab|
	def bra(bra):
		"""
		"""
		if len(bra) == 2: return (bra[0], bra[1])
		elif bra.startswith("sp"): return ("sp", bra[2:])
		else: return (bra[0], "sp")

	## Returns unique additive factors of ij
	# @param ij Number to decompose into factors
	def factor2(ij):
		"""
		Returns the unique factors for a given angular momentum value
		"""
 		return [(ij-j,j) for j in range(0, ij/2+1)]

	## Returns a list of orbitals unique with respect to permutation given an angular momentum L
	# @param L Angular momentum
	def unique_orbitals(L):
		""" 
		Returns a list of unique orbitals given angular momentum L
		"""
		o = []
		for i in range(L,-1,-1):
			for (j,k) in Shell.factor2(L-i):
				if i < j: continue
				o.append((i,j,k))
		return tuple(o) 

	## Returns unique permutations of numbers i, j, k
	# @param i First number
	# @param j Second number
	# @param k Third number
	def permute3(i,j,k):
		"""
		Returns all permutations of a given orbital
		"""
		if i == j and j == k: return [(i,j,k),]
		elif j == k: return [(i,j,k),(j,i,k),(k,j,i)]
		elif i == j: return [(i,j,k),(i,k,j),(k,j,i)]
		else: return [(i,j,k),(i,k,j),(j,i,k),(k,i,j),(j,k,i),(k,j,i)]
	
	## Returns a list of all orbitals for a given shell type
	# @param type Shell type
	def orbitals(type):
	 	"""
	 	Returns a list of orbitals for a given angular momentum.
	 	L - angular momentum
	 	"""
		o = []
		for t in type:
			L = Shell.ang_mom(t)
			for (i,j,k) in Shell.unique_orbitals(L):
				o += Shell.permute3(i,j,k) 
		return tuple(o)
	
	## Returns shell index
	# @param shell Shell type	
	def index(shell):
		if(isinstance(shell, str)):
			if(shell == "sp"):
				return Shell.index("s")
			return Shell.index(Shell.ang_mom(shell))
		L = shell
		return (L*(L + 1)*(L + 2))/6

	## Returns orbital index
	# @param (l, m, n) Oribital tuple
	def orbital_index((l,m,n)):
		L = l + m + n
		o = list(Shell.orbitals(L))
		return o.index((l,m,n))
	
	## Returns double factorial of the number n
	# @param n A number
	def factorial2(n):
		if n == 0 or n ==1:
			return 1
		else:
			return n * Shell.factorial2(n-2)
	
	## Returns orbital normalization constant
	# @param (i, j, k) Orbital tuple	 
	def orbital_norm((i,j,k)):
		n = max(i + j + k, 1)
		nx = max(i, 1)
		ny = max(j, 1)
		nz = max(k, 1)
		denom = Shell.factorial2(2*nx -1 )*Shell.factorial2(2*ny -1 )*Shell.factorial2(2*nz -1 )
		return sqrt(Shell.factorial2(2*n - 1)/denom)

	def index2(shell):
		"""
		Returns the first and last indices of a shell
		Eg. 1 and 3 for a p shell
		"""
		if(isinstance(shell,str)):
			if(shell == "sp"):
				return (Shell.index("s"), Shell.index(Shell.ang_mom("p")+1) - 1)
			return Shell.index2(Shell.ang_mom(shell))
		L = shell
		return (Shell.index(L), Shell.index(L+1) - 1)

	## Returns pairs of orbitals a and b
	# @param a First orbital
	# @param b Second orbital
	def orbitals2(a, b):
		"""
		Returns pairs of orbitals
		"""
		oa = Shell.orbitals(Shell.ang_mom(a))
		ob = Shell.orbitals(Shell.ang_mom(b))
		return tuple([(j,i) for i in ob for j in oa])

	## Finds sub-tuple t in tuple s
	# @param s Tuple
	# @param t Sub-tuple
	def find_tuple(s, t):
		"""
		"""
		for i in s:
			if i[:len(t)] == t: return i
		return None
	
	## Partitions tuple t into blocks
	# @param t Tuple to partition
	# @param block Blocking factor
	def block(t, block):
		"""
		"""
		blocks = []
		for i in range(0, len(t), block):
			blocks.append(t[i:i+block])
		return tuple(blocks)
	
	# declarations of the static method
	is_sp = staticmethod(is_sp)
	type = staticmethod(type)
	ang_mom = staticmethod(ang_mom)
	size = staticmethod(size)
	range = staticmethod(range)
	bra = staticmethod(bra)
	factor2 = staticmethod(factor2)
	unique_orbitals = staticmethod(unique_orbitals)
	permute3 = staticmethod(permute3)
	orbitals = staticmethod(orbitals)
	orbitals2 = staticmethod(orbitals2)
	index = staticmethod(index)
	index2 = staticmethod(index2)
	orbital_norm = staticmethod(orbital_norm)
	factorial2 = staticmethod(factorial2)
	find_tuple = staticmethod(find_tuple)
	block = staticmethod(block)

	def encode(seq, width):
	    e = 0
	    for (i,s) in enumerate(seq): e = e | (s << i*width)
	    return e
	encode = staticmethod(encode)

	## Creates object
	# Tuple is an immutable data type and object derived from tuple must be initialized by __new__
	# @param cls Object class
	# @param type Shell type
	def __new__(cls, type):
		"""
		"""
		return tuple.__new__(cls, Shell.orbitals(type))

	## Constructor
	# @param self Self
	# @param type Shell type
	def __init__ (self, type):
		"""
		constructor for the Shell class
		type = type of the shell
		"""
		
		## Shell type
		self.type = type
		## Shell angular momentum
		self.L = Shell.ang_mom(type)
		## sp flag
		self.sp = Shell.is_sp(type)
		## Minimum angular momentum
		self.L0 = self.L
		if self.sp: self.L0 = self.L0 - 1
		## Shell index
		self.first = Shell.index(type)
		self.hybrid = self.sp
		self.nc = 1 << self.sp

	## Returns shell as a string
	# @param self Self
	def __str__(self):
		"""
		Prints out the shell type
		"""
		return str(self.type)

	## Partitions the shell into blocks
	# If block is given, the shell is partitioned into blocks of size block. Otherwise, the shell is partitioned by orbital permutations.
	# @param self Self
	# @param block Block factor
	def partition(self, block = None):
		""" 
		Partition the orbitals. 
		If block is given the orbitals are partitioned into blocks of size block. 
		If block is None, the orbitals are partitioned by permutation
		block - blocking factor
		returns a tuple of  partitioned orbitals
		"""
		if block:
			return tuple([self[i:i+block]  for i in range(0, len(self), block)])
		else:
			blist = [[self[0]],]
			for t in self[1:]:
				if sorted(blist[-1][0]) == sorted(t): 
					blist[-1] += [t]
				else: 
					blist += [[t],]
		
		return tuple(blist)

	## Returns shell bitmask
	# @param self Self
	def type1(self): return ((-1)**self.hybrid)*self.L

	## Returns shell as an uppercase shell
	# @param self Self
	def upper(self): return str(self).upper()


## Shell pair
class Shell2(tuple):

	## Creates shell pair object
	# @param cls Object class
	# @param A First shell
	# @param B Second shell
	def __new__(cls, A, B):
		"""
		"""
		return tuple.__new__(cls, [(i,j) for j in B for i in A])

	## Initializes shell pair object
	# @param self Self
	# @param A First shell
	# @param B Second shell
	def __init__ (self, A, B):
		"""
		constructor for the Shell class
		type = type of the shell
		""" 
		## First shell
		self.A = A 
		## Second shell
		self.B = B
	
	## Returns shell pair object as a string
	# @param self Self
	def __str__(self):
		"""
		Prints out the shell type
		"""
		return str(self.A)+str(self.B)

	## Returns shell pair object as uppercase string
	# @param self Self
	def upper(self):
		return str(self).upper()

	## Partitions the shell pair orbitals into blocks
	# @param self Self
	# @param block Blocking parameter
	def partition(self, block = 1):
		""" 
		Partition the orbitals. 
		If block is given the orbitals are partitioned into blocks of size block. 
		block - blocking factor
		returns a tuple of  partitioned orbitals
		"""
		return tuple([self[i:i+block]  for i in range(0, len(self), block)])


