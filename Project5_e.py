# Project 5, part d)
# Implementation of analytic solution and makefile script for cpp code

from math import *
import numpy as np
import matplotlib.pyplot as plt

import os

os.system('g++ Project5_d.cpp -o Project5_d.o -O3 -I /Users/Henrik/FYS4150\ -\ Computational\ physics/armadillo-5.500.2/include -lblas -llapack')
os.system('./Project5_d.o Output1.txt Output2.txt Output3.txt')


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
	u = []; calc_times = []
	flag = 0
	for line in lines:
		words = line.split()
		if len(words) > 1:
			if words[0][0] == 'S':
			  T_list.append(float(words[-1]))
			  prev_line_was_T = 1
			  u_list.append(u)
			  u = []
			else: calc_times.append(float(words[-1]))
		else:
			if prev_line_was_T:
				prev_line_was_T = 0
			u.append(float(words[-1]))
	u_list.append(u)
	# First element of u_list is now empty:
	return T_list, u_list[1:], calc_times

ts_Forward, us_Forward, calc1 = read_data('Output1.txt')
ts_Backward, us_Backward, calc2 = read_data('Output2.txt')
ts_CrankNic, us_CrankNic, calc3 = read_data('Output3.txt')

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
	ax.plot(x, us_Forward[i], '--', label='$FE$')
	ax.plot(x, us_Backward[i], '--', label='$BE$')
	ax.plot(x, us_CrankNic[i], '--', label='$CN$')

#ax.legend(loc='best',fancybox='True',shadow='True')
ax.set_xlabel('Position, $x \in [0,1]$')
ax.set_ylabel('Solution, $u(x,t)$')
ax.set_title('Solution to diffusion problem')
ax.legend(loc='best',fancybox='True',shadow='True')
ax.grid()
plt.savefig('Violate_Dt.eps', format='eps', dpi=1000)
plt.show()
 
# for j in range(len(ts_Forward)):
# 	t_val = ts_Forward[j]
# 	print "Table of average absolute errors. t = %6.3f" % t_val
# 	av_er1 = 0; av_er2 = 0; av_er3 = 0
# 	for i in range(len(x)-1):
# 		x_val = x[i]
# 		exact = u_analytic(x_val, t_val)
# 		av_er1 += abs(us_Forward[j][i]-exact)
# 		av_er2 += abs(us_Backward[j][i]-exact)
# 		av_er3 += abs(us_CrankNic[j][i]-exact)
# 		#print us_Forward[j][i], exact
# 	av_er1 /= float(len(x)-1); av_er2 /= float(len(x)-1); av_er3 /= float(len(x)-1); 
# 	print " FE: %6.5f, BE: %6.5f, CN: %6.5f" % (av_er1, av_er2, av_er3)
# 	print " "



