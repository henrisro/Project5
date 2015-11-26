# Project 5, part d)
# Implementation of analytic solution and makefile script for cpp code

from math import *
import numpy as np
import matplotlib.pyplot as plt

import os

os.system('g++ Project5_d.cpp -o Project5_d.o -O3 -I /Users/Henrik/FYS4150\ -\ Computational\ physics/armadillo-5.500.2/include -lblas -llapack')
os.system('./Project5_d.o Output1_compare.txt Output2_compare.txt Output3_compare.txt')


# N is the truncation limit in the closed form expression sum:
def u_analytic(x, t, N=501):
	s = 1-x
	S = 0.0
	for n in range(1,N):
		S -= 2.0/(pi*n)*sin(n*pi*x)*exp(-n*pi*n*pi*t)
	return s + S

def read_data(filename):
	infile = open(filename, 'r')
	T_list = []
	u_list = [] # Nested list containing lists of solutions for all values of T
	lines = infile.readlines()
	prev_line_was_T = 0
	u = []
	flag = 0
	for line in lines:
		words = line.split()
		if len(words) > 1:
			T_list.append(float(words[-1]))
			prev_line_was_T = 1
			u_list.append(u)
			u = []
		else:
			if prev_line_was_T:
				prev_line_was_T = 0
			u.append(float(words[-1]))
	u_list.append(u)
	# First element of u_list is now empty:
	return T_list, u_list[1:]

ts_Forward, us_Forward = read_data('Output1_compare.txt')
ts_Backward, us_Backward = read_data('Output2_compare.txt')
ts_CrankNic, us_CrankNic = read_data('Output3_compare.txt')

N_x = len(us_Forward[0])
N_x_analytic = 501

x = np.linspace(0, 1, N_x)
x2 = np.linspace(0, 1, N_x_analytic)
u_anal = np.zeros(N_x_analytic)

plt.rcParams.update({'font.size': 14})
fig, ax = plt.subplots(1)
# Loop over time and plot solutions to see time evolution:
for i in range(len(ts_Forward)):
	t_val = ts_Forward[i]
	for j in range(N_x_analytic):
		x_val = x2[j]
		#print t_val, x_val, u_analytic(x_val, t_val)
		u_anal[j] = u_analytic(x_val, t_val)
	ax.plot(x2, u_anal, '-', label='$u(x,t=%.2f)$' % t_val) #label='t = %.2f' % t_val)
	ax.plot(x, us_Forward[i], '.--', label='$FE$')
	ax.plot(x, us_Backward[i], '.--', label='$BE$')
	ax.plot(x, us_CrankNic[i], '.--', label='$CN$')

ax.set_xlabel('Position, $x \in [0,1]$')
ax.set_ylabel('Solution, $u(x,t)$')
ax.set_title('Diffusion problem for $t_1$ and $t_2$ for three schemes.') 
#ax.plot([0.1985,0.2010],[0.525,0.525],'k-',label='$Zoom$')
#ax.plot([0.1985,0.2010],[0.530,0.530],'k-')
#ax.plot([0.1985,0.1985],[0.525,0.530],'k-')
#ax.plot([0.2010,0.2010],[0.525,0.530],'k-')

ax.plot([0.15,0.25],[0.4,0.4],'k-',label='$Zoom$')
ax.plot([0.15,0.25],[0.6,0.6],'k-')
ax.plot([0.15,0.15],[0.4,0.6],'k-')
ax.plot([0.25,0.25],[0.4,0.6],'k-')
ax.legend(loc='best',fancybox='True',shadow='True')
#ax.set_xlim(0.15,0.25)
#ax.set_ylim(0.40,0.60)
ax.grid()
plt.savefig('Comparison_n10.eps', format='eps', dpi=1000)
plt.show()

outfile = open('Results_n10.txt','w')
for j in range(len(ts_Forward)):
	t_val = ts_Forward[j]
	outfile.write("Table of average relative errors. t = %6.2f \n" % t_val)
	av_er1 = 0; av_er2 = 0; av_er3 = 0
	M = float(len(range(1,len(x)-1)))
	outfile.write("Averaging over %.1f internal points: \n" % M)
	for i in range(1,len(x)-1):
		x_val = x[i]
		#print x_val
		exact = u_analytic(x_val, t_val)
		av_er1 += abs(us_Forward[j][i]-exact)/exact
		av_er2 += abs(us_Backward[j][i]-exact)/exact
		av_er3 += abs(us_CrankNic[j][i]-exact)/exact
		#print us_Forward[j][i], exact
	av_er1 /= M; av_er2 /= M; av_er3 /= M; 
	outfile.write(" FE: %15.8f, BE: %15.8f, CN: %15.8f \n" % (av_er1, av_er2, av_er3))
	outfile.write("\n")




