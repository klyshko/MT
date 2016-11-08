__author__ = 'klyshko'

'''Program for drawing plots'''

from matplotlib import *
import numpy as np
from pylab import *
from sys import argv
import operator

pss = [50, 100, 150, 200, 230, 250, 260]
Y = []
counter = 0

for ps in pss:
	counter += 1
	folder = "{}ps/energies".format(ps)
	energies = []
	for j in range(0,10):
		energies.append("{}{}.dat".format(folder, j))
	dic = {}
	for f in energies:
		fi = open(f)
		for line in fi:
			words = line.split()
			while words[1] != 'nan':
				dic[f] = int(words[0]) #* ps / 10.0**12
				break
	sorted_dic = sorted(dic.items(), key=operator.itemgetter(1))
	print sorted_dic
	y = []
	for item in sorted_dic:
		y.append(item[1])
	Y.append(y)


for y in Y:
	y = np.array(y)	

X = np.arange(0.1, 1.1, 0.1)

output = open("probabilities.dat", "w")
output.write("prob      \t50ps      \t100ps      \t150ps     \t200ps      \t230ps      \t250ps\t    260ps\t\n")

for i in range(0, len(X)):
	output.write("%f" % X[i])
	for k in range(0,len(pss)):
		output.write("\t%9d" % Y[k][i])
	output.write("\n")

output.close()
print("file saved")

"""
figure(figsize=(18, 12))
xlabel("time (s)")
ylabel("pobability")
xlim(0, 0.1)
title("timestep = {} ps".format(ps))
for k in range(0,len(pss)):
	plot(Y[k], X, "-", linewidth=2)
legend(['50ps', '100ps', '150ps', '200ps', '230ps', '250ps', '260ps'], loc='upper right')
#plot(Y[0], X, "-", Y[1], X, "-", Y[2], X, "-", Y[3], X, "-", Y[4], X, "-", Y[5], X, "-", Y[6], X, "-")
print "plot saved to prob.png"
savefig("prob.png")
grid()
show()
"""



