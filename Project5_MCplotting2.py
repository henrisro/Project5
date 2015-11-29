# Project 5, part g), plotting part. 
# This program is only used for plotting multiple histograms in the
# same window. Figures produced here is used in the report.

import os
import matplotlib.pyplot as plt
import numpy as np
from math import pi, sqrt

# First pair of plots:
name1 = 'MC_gaussian_lowT1.txt'
name2 = 'MC_gaussian_highT1.txt'

# Final pair of plots:
name3 = 'MC_gaussian_lowT2.txt'
name4 = 'MC_gaussian_highT2.txt'

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

tsteps1, walkers1, dt1, N_final1, calctime1, positions1, varlist1 = read_data_MC(name3)
tsteps2, walkers2, dt2, N_final2, calctime2, positions2, varlist2 = read_data_MC(name4)

low_T = dt1*tsteps1
high_T = dt2*tsteps2

l01 = sqrt(2.0*dt1)
m = int(1/l01)

# Calculating the analytic solution:
N_x_analytic = 500
x_anal = np.linspace(0, 1, N_x_analytic)
u_anal_low = np.zeros(N_x_analytic); u_anal_high = np.zeros(N_x_analytic)
for j in range(N_x_analytic):
	x_val = x_anal[j]
	u_anal_low[j] = u_analytic(x_val, low_T)
	u_anal_high[j] = u_analytic(x_val, high_T)

bins = [i*l01 for i in range(m+1)]

fig, ax = plt.subplots(1)
plt.rcParams.update({'font.size': 14})
n1, bins, patches1 = plt.hist(positions1, bins, normed=1) #, facecolor='red', alpha=0.75)
plt.clf()
n2, bins, patches2 = plt.hist(positions2, bins, normed=1) #, facecolor='red', alpha=0.75)
plt.clf()

new_n1 = []; new_n2 = []
norm1 = float(walkers1)/(walkers1+N_final1)/2.0
norm2 = float(walkers2)/(walkers2+N_final2)/2.0
for i in range(len(n1)):
	val1 = n1[i]; val2 = n2[i]
	xval = bins[i]+l01/2.0
	yval1 = -float(val1)*norm1 + 1.0 - xval
	yval2 = -float(val2)*norm2 + 1.0 - xval
	new_n1.append(yval1); new_n2.append(yval2)

#Adding u_s(x): the stationary solution.
for j in range(N_x_analytic):
	x_val = x_anal[j]
	u_anal_low[j] += 1 - x_val
	u_anal_high[j] += 1 - x_val

width = bins[1] - bins[0] # Width of each bin.
x1 = np.ravel(zip(bins[:-1], bins[:-1]+width))
y1 = np.ravel(zip(new_n1, new_n1))
y2 = np.ravel(zip(new_n2, new_n2))

plt.plot(x_anal, u_anal_low, 'k-', label='$u(x,t=%.2f$)' % low_T)
plt.plot(x1, y1,'r-',label='MC, $N_{MC} = %.0E$' % walkers1)
plt.plot(x_anal, u_anal_high, 'b-', label='$u(x,t=%.2f$)' % high_T)
plt.plot(x1, y2,'g-',label='MC, $N_{MC} = %.0E$' % walkers2)
plt.xlabel('Position, $x\in [0,1]$')
plt.ylabel('Solution, $u(x,t)$')
plt.legend(loc='best',fancybox='True',shadow='True')
plt.title('Monte Carlo simulation of diffusion problem. \n Normal distributed step length.')
plt.grid('on')
plt.savefig('MC_gaussian_twoplot_T2.eps', format='eps', dpi=1000)
plt.show()

