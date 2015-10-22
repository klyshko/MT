import matplotlib
from pylab import plot, show, ylim, yticks, grid, legend, title, xlabel, ylabel, savefig, subplot
from numpy import sin, cos, exp, pi, arange
from sys import argv

def pot(D, A, scale, sigma, x):
	return (D*(1-exp(-A*(x-4)))**2 - D + scale*(sigma/x)**6)

def dpot(D, A, scale, sigma, x):
	return -( 2*A*D*(1-exp(-A*(x-4)))*exp(-A*(x-4)) + scale*(sigma**6)/(x)**7)

def read_config(filename):
    config = {}
    for i in open(filename, 'r').readlines():
        if i[0] != '#' and i[0] != '\n' and i[0] != ' ':
            t = i.split(None, 1)
            config[t[0]] = t[1].replace('\n','').split('#')[0]
    return config

def config_replace(params, x):
    for par in params: x = x.replace('<'+par+'>', params[par])
    return x

def config_get_replaced(params, a):
    return config_replace(params, params[a])

p_ = read_config(argv[1])

A_long = float(config_get_replaced(p_, "A_long"))
D_long = float(config_get_replaced(p_, "D_long"))
A_lat = float(config_get_replaced(p_, "A_lat"))
D_lat = float(config_get_replaced(p_, "D_lat"))

lj_scale = float(config_get_replaced(p_, "LJScale"))
lj_sigma = float(config_get_replaced(p_, "LJSigma"))


t = arange(3.93, 5.9, 0.01)

matplotlib.pyplot.figure(num=None, figsize=(16, 12), dpi=300, facecolor='w', edgecolor='k')
subplot(221)
title('Longitudinal')
s1 = pot(D_long, A_long, lj_scale, lj_sigma, t)
#ylim(-20,10)
plot(t, s1)
grid()
xlabel('distance, nm')
ylabel('energy, kcal/mol')

subplot(222)
title('Lateral')
s1 = pot(D_lat, A_lat, lj_scale, lj_sigma, t)
#ylim(-20,10)
plot(t, s1)
grid()
xlabel('distance, nm')
ylabel('energy, kcal/mol')

subplot(223)
s1 = dpot(D_long, A_long, lj_scale, lj_sigma, t)
plot(t, s1)
grid()
xlabel('distance, nm')
ylabel('force, kcal/(mol*nm)')

subplot(224)
s1 = dpot(D_lat, A_lat, lj_scale, lj_sigma, t)
plot(t, s1)
grid()
xlabel('distance, nm')
ylabel('force, kcal/(mol*nm)')

savefig("plot_morse_lj.png", dpi=300, format='png')
#show()
