import matplotlib
from pylab import plot, show, ylim, yticks, grid, legend, title, xlabel, ylabel, savefig, subplot
from numpy import sin, cos, exp, pi, arange
from sys import argv

def pot(D, A, a, r, w, scale, sigma, x):
	return (D*(1-exp(-A*(x-4)))**2 - D + a*exp(-((x-4)-r)**2/(w*w)) + scale*(sigma/x)**6)

def dpot(D, A, a, r, w, scale, sigma, x):
	return -( 2*A*D*(1-exp(-A*(x-4)))*exp(-A*(x-4)) - a*exp(-((x-4)-r)**2/(2*w*w))*((x-4)-r)/(w*w) + scale*(sigma**6)/(x)**7)

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

a_lat_b = float(config_get_replaced(p_, "a_barr_lat"))
r_lat_b = float(config_get_replaced(p_, "r_barr_lat"))
w_lat_b = float(config_get_replaced(p_, "w_barr_lat"))

a_long_b = float(config_get_replaced(p_, "a_barr_long"))
r_long_b = float(config_get_replaced(p_, "r_barr_long"))
w_long_b = float(config_get_replaced(p_, "w_barr_long"))

lj_scale = float(config_get_replaced(p_, "LJScale"))
lj_sigma = float(config_get_replaced(p_, "LJSigma"))


t = arange(3.75, 6.0, 0.01)

matplotlib.pyplot.figure(num=None, figsize=(16, 12), dpi=300, facecolor='w', edgecolor='k')
subplot(221)
title('Longitudinal')
s1 = pot(D_long, A_long, a_long_b, r_long_b, w_long_b, lj_scale, lj_sigma, t)
#addcode
longfile = open("long.dat","w")
for tt in t:
    longfile.write("{}  {}\n".format(tt, pot(D_long, A_long, a_long_b, r_long_b, w_long_b, lj_scale, lj_sigma, tt)))
longfile.close()
#addcode
plot(t, s1)
grid()
xlabel('distance, nm')
ylabel('energy, kcal/mol')

subplot(222)
title('Lateral')
s1 = pot(D_lat, A_lat, a_lat_b, r_lat_b, w_lat_b, lj_scale, lj_sigma, t)
#ylim(-20,10)
#addcode
latfile = open("lat.dat","w")
for tt in t:
    latfile.write("{}  {}\n".format(tt, pot(D_lat, A_lat, a_lat_b, r_lat_b, w_lat_b, lj_scale, lj_sigma, tt)))
latfile.close()
#addcode
plot(t, s1)
grid()
xlabel('distance, nm')
ylabel('energy, kcal/mol')

subplot(223)
s1 = dpot(D_long, A_long, a_long_b, r_long_b, w_long_b, lj_scale, lj_sigma, t)
plot(t, s1)
grid()
xlabel('distance, nm')
ylabel('force, kcal/(mol*nm)')

subplot(224)
s1 = dpot(D_lat, A_lat, a_lat_b, r_lat_b, w_lat_b, lj_scale, lj_sigma, t)
plot(t, s1)
grid()
xlabel('distance, nm')
ylabel('force, kcal/(mol*nm)')

savefig("plot_morse_bar_lj.png", dpi=300, format='png')
#show()
