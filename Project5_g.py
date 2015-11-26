# Project 5, part g): 
# Implementation of random walk as diffusion problem
# This is a makefile for Project5_g.cpp

import os
import matplotlib.pyplot as plt
import numpy as np
from math import pi, sqrt

name = 'MC_gaussian_highT2.txt'
outfilename = 'Benchmark_highT2_gaussian.txt'

os.system('g++ -c Project5_g.cpp lib.cpp')
os.system('g++ -o Project5_g.out Project5_g.o lib.o -O3')
# Output file should be provided here as well:
os.system('./Project5_g.out %s 100000' % name)

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
	dt = float(lines[2].split()[-1])

	for line in lines[3:-2]:
		words = line.split()
		position_list.append(float(words[0]))
		var_list.append(float(words[1]))

	N_final = int(lines[-2].split()[-1])
	calctime = float(lines[-1].split()[-1])
	return time_steps, walks, dt, N_final, calctime, position_list, var_list

tsteps, walkers, dt, N_final, calctime, positions, varlist = read_data_MC(name)

# Calculating standard deviation and using it as an estimate of the error:
average_var = 0.0
for k in range(walkers):
 	average_var += varlist[k]
std_N = sqrt(average_var)/sqrt(walkers+N_final)/sqrt(tsteps)
print " Average std square: ", std_N

l0 = sqrt(2.0*dt)
m = int(1/l0)
print "tsteps: ", tsteps, " walkers: ", walkers, " l0: ", l0, " and time step:", dt
time = dt*tsteps
print "Final time: ", time
bins = [i*l0 for i in range(m+1)]

# Calculating the analytic solution:
N_x_analytic = 500
x_anal = np.linspace(0, 1, N_x_analytic)
u_anal = np.zeros(N_x_analytic)
for j in range(N_x_analytic):
	x_val = x_anal[j]
	u_anal[j] = u_analytic(x_val, time)

# Creating histogram and normalizing it:
fig, ax = plt.subplots(1)
plt.rcParams.update({'font.size': 14})
n, bins, patches = plt.hist(positions, bins, normed='True') #, facecolor='red', alpha=0.75)
plt.clf()
new_n = []
norm2 = 1.0/(max(n)/abs(min(u_anal))) # This is only for visual comparison.
# Correct norming:
norm = float(walkers)/(walkers+N_final)/2.0 # Correct normalization
print " Norm = ", norm
print " Ratio between amplitudes for comarison: ", norm2
for i in range(len(n)):
	val = n[i]
	xval = bins[i]+l0/2
	yval = -float(val)*norm + 1.0 - xval
	new_n.append(yval)

# Adding u_s(x): the stationary solution for visualization:
for j in range(N_x_analytic):
	x_val = x_anal[j]
	u_anal[j] += 1 - x_val

# Making beautiful plots:
width = bins[1] - bins[0] # Width of each bin.
x = np.ravel(zip(bins[:-1], bins[:-1]+width))
y = np.ravel(zip(new_n, new_n))
plt.plot(x_anal, u_anal, 'k-', label='$u(x,t=%.3f$)' % time)
plt.plot(x, y,'b-',label='Monte Carlo $N_{MC} = %d$' % walkers)
plt.xlabel('Position, $x\in [0,1]$')
plt.ylabel('Solution, $u(x,t)$')
plt.legend(loc='best',fancybox='True',shadow='True')
plt.title('Monte Carlo simulation of diffusion problem. \n Normal distributed step length.')
plt.grid('on')
plt.savefig(name[:-4]+'.eps', format='eps', dpi=1000)
plt.show()

# Writing the average (absolute) error to a txt file:
outfile = open(outfilename,'w')
outfile.write("Table of average relative errors. t = %6.2f, N_MC = %d, dt = %8.5f, tsteps = %d, average std. = %10.8f \n" % (time, walkers, dt, tsteps, std_N))
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