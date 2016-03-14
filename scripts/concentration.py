"""10 steps to run assembly simulations with constant concentration:

1) If the necessary concentration is C [microMole], then it equals 6*C*1e-7 [1/nm^3]. Put C to cond.conf for final simulations
2) Find your Volume parameters: rep_h and rep_r in "morse.conf". 
3) Calculate volume V = 3.14 * rep_r * rep_r * rep_h  (usually r = 80 nm, h = 160 nm)
4) Calculate necessary # of free particles: # = 6*C*1e-7 * V
5) Generate MT with make_mt_cncnt.py file xyz_#_?.pdb and ang_#_?.pdb, where '?' is any number > #
6) Compile with  LJ and REPULSIVE, no MORSE, no ASSEMBLY, no MT_LENGTH, no CONCentration. 
7) Simulate *.pdb files for 10 million steps and fix = 2, runnum = 1;
8) Take result.pdb as initial structure
9) Compile with LJ, CONC, MORSE, ASSEMBLY, MT_LENGTH
10 Simulate trajectory with any runnum, fix = 2 and so on.

11) PROFIT!!

"""

from sys import argv
from math import cos,sin,pi

if len(argv) != 4:
		print 'Usage: %s concentration (in [muMol/L]) rep_r rep_h \n'% argv[0]
		conc  = float(argv[1])
		rep_r = 80.0
		rep_h = 160.0

		V = pi * rep_r **2 * rep_h
		C = 6 * conc * 1e-7
		print "conc = %f [muMol/L] = %f [dimers/nm^3]\n" % (conc,C)
		print "Volume = %f [nm^3]  - used r = 80, h = 160 by default\n " % V
		print "Free particles # (dimers) = %d \n" % int(C * V) 
		print "Launch python make_mt_cncntr.py %d extra#\n" % ( 2 * (int(C * V / 13) + 1))
else:		
		conc  = float(argv[1])
		rep_r = float(argv[2])
		rep_h = float(argv[3])

		V = pi * rep_r **2 * rep_h
		C = 6 * conc * 1e-7
		print "\nconc = %f [muMol/L] = %f [dimers/nm^3]\n" % (conc,C)
		print "Volume = %f [nm^3]\n " % V
		print "Free particles # (dimers) = %d \n" % int(C * V) 
		print "Launch python make_mt_cncntr.py %d extra#\n" % ( 2 * (int(C * V / 13) + 1))
		
print """6) Compile with  LJ and REPULSIVE, no MORSE, no ASSEMBLY, no MT_LENGTH, no CONCentration. 
7) Simulate *.pdb files for 10 million steps and fix = 2, runnum = 1;
8) Take result.pdb as initial structure	
9) Compile with LJ, CONC, MORSE, ASSEMBLY, MT_LENGTH
10) Simulate trajectory with any runnum, fix = 2 and so on."""