from pdb import PDBAtom
from sys import argv
from math import cos,sin,pi,atan
r_mon      = 2.0
R_MT   	   = 8.12
alpha_step = 2*pi/13.0
z_step	   = 6.0/13.0*r_mon
tail_step  = r_mon/2
theta_tail = 0.21
if len(argv) != 3 and len(argv) != 4:
		print 'Usage: %s MT_length tail_length [-d]\nAll lengths should be specified in monomers.\nUse -d to delete "tail" of one protofilament.'%argv[0]
		exit(0)
MT_len   = int(argv[1])
tail_len = int(argv[2])

if len(argv) == 4 and argv[3] == '-d':
		delete = True
else:
		delete = False
if MT_len<tail_len:
		print 'Tail can not be longer, than MT'
		exit(0)

pdbxyzfilename = 'structs/xyz_{}.pdb'.format(MT_len)
pdbangfilename = 'structs/ang_{}.pdb'.format(MT_len)
pdbxyz = open(pdbxyzfilename, "w")
pdbang = open(pdbangfilename, "w")

pdbstr1   = 'ATOM     14   CB ALA A   7       8.120  -0.000  52.923  0.00  8.12         A'
pdbat1    = PDBAtom(pdbstr1)
pdbstr2   = 'ATOM     14   CB GLY A   7       8.120  -0.000  52.923  0.00  8.12         A'
pdbat2    = PDBAtom(pdbstr2)
chains   = [chr(i) for i in range(ord('A'), ord('M')+1)]
for chain in chains:
		for mon in range(MT_len):	
			if (MT_len - mon) < tail_len:
					if delete == True and chain == chains[len(chains) - 1]:
							break
					R 		= R_MT + (tail_len - (MT_len - mon))*2*r_mon*sin(theta_tail) #tail_step*(tail_len - (MT_len - mon))
					pdbat2.z = theta_tail #atan(tail_step/z_step)
					pdbat1.z = 2*r_mon*mon + z_step*(ord(chain) - ord('A'))
			else:
					R 		 = R_MT
					pdbat2.z  = 0.0
					pdbat1.z  = 2*r_mon*mon + z_step*(ord(chain) - ord('A'))
			alpha = alpha_step*(ord(chain) - ord('A'))
			pdbat2.y = -alpha
			pdbat2.x = 0.0
			x = R*cos(-alpha)
			y = R*sin(-alpha)
			pdbat1.chainid    = chain
			pdbat2.chainid    = chain
			pdbat1.x          = x
			pdbat1.y          = y
			pdbat1.residueseq = mon/2 + 1
			pdbat2.residueseq = mon/2 + 1
			pdbat1.number     = mon + MT_len*(ord(chain)-ord(chains[0])) + 1 
			pdbat2.number     = mon + MT_len*(ord(chain)-ord(chains[0])) + 1 
			if mon%2 == 0:
					pdbat1.name = 'CA'
					pdbat2.name = 'CA'
			else:
					pdbat1.name = 'CB'
					pdbat2.name = 'CB'
			pdbxyz.write("%s\n" % str(pdbat1))
			pdbang.write("%s\n" % str(pdbat2))
pdbxyz.write('END')
pdbang.write('END')
pdbxyz.close()
pdbang.close()
