__author__ = 'klyshko'

'''Program for drawing plot of energies'''

from matplotlib import *
import numpy as np
from pylab import *
from sys import argv

filename = argv[1]


f = open(filename,'r')
x = []
y = [[] for i in range(7)]

for line in f:
	if not (line[0] == '@' or line[0] == '#' or line[0] == '&'):
   		z = line.split()
   		x.append(float(z[0]))
   		for i in range(1,8):
   			y[i-1].append(float(z[i]))

x = np.array(x)

xmin = min(x)
xmax = max(x)
X = np.linspace(xmin, xmax, num=10000)

harm = np.array(y[0])
longit = np.array(y[1])
lat = np.array(y[2])
psi = np.array(y[3])
phi = np.array(y[4])
theta = np.array(y[5])
lj = np.array(y[6])

plot(x, harm, "y-", x, longit, "g.", x, lat, "b.", x, theta, "r.", x, lj, "g-")
newfile = "%s.png" % filename[:-4]
print "plot saved to", newfile
savefig(newfile)
grid()
show()

"""
#plot2 = figure(2)
scatter(r,Q, s=0.1)
plot(X, np.array(map(F,X)), "-", color = "g")
axis([0, 0.31, 0, 1.2])
savefig("3.png")
show()

def poly(x, b,c,d):
    return d*(x+b)**c
"""

#raw_input()
