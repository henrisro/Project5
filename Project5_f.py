# Project 5, part f): 
# Implementation of random walk as diffusion problem.
# Makefile for Project5_f.cpp

import os
import matplotlib.pyplot as plt
import numpy as np
from math import pi, sqrt

name = 'MC_uniform_highT2.txt'
outfilename = 'Benchmark_highT2.txt'

os.system('g++ -c Project5_f.cpp lib.cpp')
os.system('g++ -o Project5_f.out Project5_f.o lib.o -O3')
# Output file should be provided here as well:
os.system('./Project5_f.out %s 100000' % name)

# N is the truncation limit in the closed form expression sum:
def u_analytic(x, t, N=501):
	S = 0.0
	for n in range(1,N):
		S -= 2.0/(pi*n)*np.sin(n*pi*x)*np.exp(-n*pi*n*pi*t)
	return S # Imlementing function as v(x) and adding u_s(x) in the end

def read_data_MC(filename):
	infile = open(filename, 'r')
	position_list = []; var_list = []

	lines = infile.readlines()
	time_steps = int(lines[0].split()[-1])
	walks = int(lines[1].split()[-1])
	l0 = float(lines[2].split()[-1])
	dt = float(lines[3].split()[-1])

	for line in lines[4:-2]:
		words = line.split()
		position_list.append(float(words[0]))
		var_list.append(float(words[1]))

	N_final = int(lines[-2].split()[-1])
	calctime = float(lines[-1].split()[-1])
	return time_steps, walks, l0, dt, N_final, calctime, position_list, var_list

tsteps, walkers, l0, dt, N_final, calctime, positions, varlist = read_data_MC(name)

average_var = 0.0
for k in range(walkers):
	average_var += varlist[k]
average_std = sqrt(average_var)/sqrt(walkers+N_final)/sqrt(tsteps)
print " Average std: ", average_std

m = int(1/l0)
print "tsteps: ", tsteps, " walkers: ", walkers, " l0: ", l0, " and time step:", dt
time = dt*tsteps
print "Final time: ", time
bins = [i*l0+l0/2.0 for i in range(m)]

# Calculating the analytic solution:
N_x_analytic = 500
x_anal = np.linspace(0, 1, N_x_analytic)
u_anal = np.zeros(N_x_analytic)
for j in range(N_x_analytic):
	x_val = x_anal[j]
	u_anal[j] = u_analytic(x_val, time)

fig, ax = plt.subplots(1)
plt.rcParams.update({'font.size': 14})
n, bins, patches = plt.hist(positions, bins, normed=1) #, facecolor='red', alpha=0.75)
plt.clf()
area_histogram = l0*sum(n)
print " Histogram area: ", area_histogram
new_n = []
norm2 = 1.0/(max(n)/abs(min(u_anal)))
# Correct norming:
norm3 = float(walkers)/(walkers+N_final)/2.0
print " Norm3 = ", norm3
print " Fraction of amplitudes: ", norm2
for i in range(len(n)):
	val = n[i]
	xval = bins[i]+l0/2.0
	yval = -float(val)*norm3 + 1.0 - xval
	new_n.append(yval)

#Adding u_s(x): the stationary solution.
for j in range(N_x_analytic):
	x_val = x_anal[j]
	u_anal[j] += 1 - x_val

width = bins[1] - bins[0] # Width of each bin.
x = np.ravel(zip(bins[:-1], bins[:-1]+width))
y = np.ravel(zip(new_n, new_n))
plt.plot(x_anal, u_anal, 'k-', label='$u(x,t=%.3f$)' % time)
plt.plot(x, y,'b-',label='Monte Carlo $N_{MC} = %d$' % walkers)
plt.xlabel('Position, $x\in [0,1]$')
plt.ylabel('Solution, $u(x,t)$')
plt.legend(loc='best',fancybox='True')
plt.title('Monte Carlo simulation of diffusion problem. \n Uniform step length.')
plt.grid('on')
plt.savefig(name[:-4]+'.eps', format='eps', dpi=1000)
plt.show()

outfile = open(outfilename,'w')
outfile.write("Table of average relative errors. t = %6.2f, N_MC = %d, dt = %8.5f, tsteps = %d, average std. = %10.8f \n" % (time, walkers, dt, tsteps, average_std))
outfile.write("Calculation time: %8.5f s \n" % calctime)
av_er = 0.0
M = float(len(new_n))
outfile.write("Averaging over %.1f internal points: \n" % M)
for i in range(len(new_n)):
	x_val = bins[i]+l0/2.0
	exact = u_analytic(x_val, time) + 1.0 - x_val
	av_er += abs(new_n[i] - exact)
	#print " Average error for x = ", x_val, " numeric: ", new_n[i], " analytic: ", exact, " abs(error): ", abs(new_n[i] - exact)
av_er /= M
outfile.write(' %10.8f ' % av_er)