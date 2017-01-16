from pdb import PDBAtom
from sys import argv
from math import cos,sin,pi,atan
r_mon      = 2.0
R_MT   	   = 8.12
alpha_step = 2*pi/13.0
z_step	   = 6.0/13.0*r_mon
tail_step  = r_mon/2
if len(argv) != 3:
		print 'Usage: %s MT_length MT_extra\nAll lengths should be specified in monomers. MT_length should be even: second half of tubule will bw for extra monomers during assembly\n'%argv[0]
		exit(0)
MT_len   = int(argv[1])
MT_extra   = int(argv[2])

if MT_len % 2:
	MT_len += 1

pdbxyzfilename = 'structs/xyz_{}_{}.pdb'.format(MT_len, MT_extra)
pdbangfilename = 'structs/ang_{}_{}.pdb'.format(MT_len, MT_extra)
pdbxyz = open(pdbxyzfilename, "w")
pdbang = open(pdbangfilename, "w")

pdbstr1   = 'ATOM     14   CB ALA A   7       8.120  -0.000  52.923  0.00  8.12         A'
pdbat1    = PDBAtom(pdbstr1)
pdbstr2   = 'ATOM     14   CB GLY A   7       8.120  -0.000  52.923  0.00  8.12         A'
pdbat2    = PDBAtom(pdbstr2)

chains   = [chr(i) for i in range(ord('A'), ord('M')+1)]
for chain in chains:
		for mon in range(MT_len):	
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

for chain in chains:
		for mon in range(MT_len, MT_len + MT_extra):	
			R 		 = R_MT
			pdbat2.z  = 0.0
			pdbat2.y = 0.0
			pdbat2.x = 0.0
			pdbat1.chainid    = 'X'
			pdbat2.chainid    = 'X'
			pdbat1.x          = 0.0
			pdbat1.y          = 0.0
			pdbat1.residueseq = mon/2 + 1
			pdbat2.residueseq = mon/2 + 1
			pdbat1.number     = MT_len * (len(chains) -1) + mon + MT_extra*(ord(chain)-ord(chains[0])) + 1 
			pdbat2.number     = MT_len * (len(chains) -1) + mon + MT_extra*(ord(chain)-ord(chains[0])) + 1 
			if mon%2 == 0:
					pdbat1.name = 'CA'
					pdbat2.name = 'CA'
					pdbat1.z  = 600.0
			else:
					pdbat1.name = 'CB'
					pdbat2.name = 'CB'
					pdbat1.z  = 604.0
			pdbxyz.write("%s\n" % str(pdbat1))
			pdbang.write("%s\n" % str(pdbat2))

pdbxyz.write('END')
pdbang.write('END')
pdbxyz.close()
pdbang.close()
