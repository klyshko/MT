__author__ = 'klyshko'

'''Program for drawing plots'''

from matplotlib import *
import numpy as np
from pylab import *
from sys import argv
import operator

ps = 50.0


filename = "energies"
energies = []
for i in range(0,10):
	energies.append("{}{}.dat".format(filename, i))

dic = {}
arr = []
for f in energies:
	fi = open(f)
	for line in fi:
		words = line.split()
		while words[1] != 'nan':
			dic[f] = float(words[0]) * ps / 10.0**12
			break
Y = []
sorted_dic = sorted(dic.items(), key=operator.itemgetter(1))

print sorted_dic

for item in sorted_dic:
	Y.append(item[1])
Y = np.array(Y)	
X = np.arange(0.1, 1.1, 0.1)


figure(figsize=(18, 12))
xlabel("time (s)")
ylabel("pobability")
title("timestep = {} ps".format(ps))
plot(Y, X, "b-")
print "plot saved to prob.png"
savefig("prob.png")
grid()
show()



