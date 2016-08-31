from pdb import *
from sys import argv
from math import cos,sin,pi,atan,sqrt

if len(argv) != 2:
		print('Usage: %s filename' % argv[0])
		exit(0)
pdbfile  = str(argv[1])

r_mon = 2.0


#atoms = PDBSystem(pdbfile)
#atoms1 = PDBReadFile(pdbfile)
trajectory = PDBTrajectory(pdbfile)
print('# of dimers / # of oligos / # of atoms /')
for frame in trajectory.frames:
	atoms = frame.atoms
	oligomers = []
	for atom1 in atoms: 
		oligomer = [atom1.number]
		for atom2 in atoms:
			dist2 = (atom1.x - atom2.x)**2 + (atom1.y - atom2.y)**2 + (atom1.z - atom2.z)**2
			if dist2 < (2.5 * r_mon)**2 and atom1 != atom2:
				#print sqrt(dist2)
				oligomer.append(atom2.number)
		oligomers.append(oligomer)

	for index1 in range(0,len(oligomers)):
		oligomer1 = oligomers[index1]
		for index2 in range(0,len(oligomers)):
			oligomer2 = oligomers[index2]
			if set(oligomer1).intersection(set(oligomer2)):
				oligomers[index1] = list(set(oligomer1).union(set(oligomer2)))
				oligomers[index2] = list(set(oligomer1).union(set(oligomer2)))

	array = [set(i) for i in oligomers]
	new = []
	for i in array:
		if i not in new:
			new.append(i)

	dimers = 0
	for i in new:
		i = sorted(list(i))
		'''if len(i) % 2:
			print 'YOBA!'
		print i'''
		if len(i) == 2:
			dimers += 1
	print('{}\t{}\t{}'.format(dimers, len(new), len(atoms)))

