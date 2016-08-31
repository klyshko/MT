from pdb import PDBAtom
from sys import argv
from math import cos,sin,pi,atan
import random

r_mon      = 2.0
rep_r = 150.
rep_h = 300.
conc = 50.

V = pi * rep_r **2 * rep_h
C = 6 * conc * 1e-7 #[dimers/nm^3]
'''

'''
n = int(V*C + 1)
print n, V, C
atoms = []

def generate(n):
	global atoms
	while len(atoms) < n:
		z = random.uniform(0.0, rep_h - 2 *r_mon)
		while True:
			x = random.uniform(-rep_r, rep_r)
			y = random.uniform(-rep_r, rep_r)
		#	for xyz in atoms:
			if x**2 + y**2 < rep_r**2:
					#and (x - xyz[0])**2 + (y - xyz[1])**2 + (z - xyz[2])**2 > r_mon*r_mon and (x - xyz[0])**2 + (y - xyz[1])**2 + (z + 2* r_mon - xyz[2])**2 > r_mon*r_mon:
				break
		atoms.append((x,y,z))
	#for atom in atoms:
		#print atom, atoms.index(atom)	

def sparse():
	global atoms		
	for atom1 in atoms:
		for atom2 in atoms:
			if atom1 != atom2 and (atom1[0] - atom2[0])**2 + (atom1[1] - atom2[1])**2 + (atom1[2] - atom2[2])**2 < 16*r_mon**2:
				atoms.remove(atom2)			

while True:
	generate(n)
	sparse()
	if len(atoms) == n:
		break

pdbxyzfilename = 'volxyz_{}.pdb'.format(int(conc))
pdbangfilename = 'volang_{}.pdb'.format(int(conc))

pdbxyz = open(pdbxyzfilename, "w")
pdbang = open(pdbangfilename, "w")

pdbstr1   = 'ATOM     14   CB ALA A   7       8.120  -0.000  52.923  0.00  8.12         A'
pdbat1    = PDBAtom(pdbstr1)
pdbstr2   = 'ATOM     14   CB GLY A   7       8.120  -0.000  52.923  0.00  8.12         A'
pdbat2    = PDBAtom(pdbstr2)

for atom in atoms:
	pdbat1.x = atom[0]
	pdbat1.y = atom[1]
	pdbat1.z = atom[2]
	pdbat1.chainid    = 'A'
	pdbat1.residueseq = atoms.index(atom) + 1
	pdbat1.number     = 2 * atoms.index(atom) + 1 
	pdbat1.name = 'CA'
	pdbxyz.write("%s\n" % str(pdbat1))

	pdbat1.x = atom[0]
	pdbat1.y = atom[1]
	pdbat1.z = atom[2]+2*r_mon
	pdbat1.chainid    = 'A'
	pdbat1.residueseq = atoms.index(atom) + 1
	pdbat1.number     = 2 * atoms.index(atom) + 2 
	pdbat1.name = 'CB'
	pdbxyz.write("%s\n" % str(pdbat1))

	pdbat2.x = 0.0
	pdbat2.y = 0.0
	pdbat2.z = 0.0
	pdbat2.chainid    = 'A'
	pdbat2.residueseq = atoms.index(atom) + 1
	pdbat2.number     = 2 * atoms.index(atom) + 1 
	pdbat2.name = 'CA'
	pdbang.write("%s\n" % str(pdbat2))

	pdbat2.x = 0.0
	pdbat2.y = 0.0
	pdbat2.z = 0.0
	pdbat2.chainid    = 'A'
	pdbat2.residueseq = atoms.index(atom) + 1
	pdbat2.number     = 2 * atoms.index(atom) + 2 
	pdbat2.name = 'CB'
	pdbang.write("%s\n" % str(pdbat2))
pdbxyz.write('END')
pdbang.write('END')
pdbxyz.close()
pdbang.close()

