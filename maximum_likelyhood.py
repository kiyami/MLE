import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import chisquare
import scipy.stats as ss
from pylab import *
import random

#########################################################################
# Asymmetric error calculation with maximum likelihood estimator (MLE)
# A and B values with asymmetric errors
# (A_neg -> minus error, A_poz -> plus error)
# Operation -> select the operation between A and B (A+B, A-B...)
# N -> number of produced random value
# delta_L -> confidence interval
# (0.5 ML interval is means 1 sigma confidence interval)
#########################################################################

N = 10000

A,A_neg,A_poz = 10.0,2.0,3.0
B,B_neg,B_poz = 5.0,0.5,0.6

graphp_limits_A = [A - 2.0*A_neg,A + 2.0*A_poz]
graphp_limits_B = [B - 2.0*B_neg,B + 2.0*B_poz]
graphp_limits_sum = [10,20]

#confidence interval
delta_L = 0.5

# select an operation (+,-,*,/)
operation = "+"

#########################################################################

def func(data,mu,sigma_neg,sigma_poz):
	a = (mu - data)
	b = (2.0*sigma_poz*sigma_neg)/(sigma_poz+sigma_neg)
	c = (sigma_poz-sigma_neg)/(sigma_poz+sigma_neg)
	e = b + (c*(-a))
	f = a/e
	chi = -0.5*((a/e)**2.0)
	return chi

#########################################################################

def mle(data,sigma_neg,sigma_poz,N):
	chi_list = []
	x_plot = np.linspace(data - 2.0*sigma_neg, data + 2.0*sigma_poz, N)
	for i in range(N):
		a = (data - x_plot[i])
		b = (2.0*sigma_poz*sigma_neg)/(sigma_poz+sigma_neg)
		c = (sigma_poz-sigma_neg)/(sigma_poz+sigma_neg)
		e = b + (c*(-a))
		f = a/e
		chi_list = np.append(chi_list,f**2.0)
	return -0.5*chi_list

#########################################################################

def find_sigma(chi_list,x_plot,N,delta_L):
	max_value = max(chi_list)
	sigma_value = max_value - delta_L
	print("max_value",max_value)
	print("1_sigma_value",sigma_value)

	for i in range(N):
		if chi_list[i] == max_value:
			max_x = x_plot[i]
	print("max_x",max_x)

	left_sigma = 0.0
	right_sigma = 0.0

	indicator = 1
	for i in range(N):
		if indicator == 1:
			if chi_list[i] >= (max_value - delta_L):
				left_sigma = x_plot[i]
				indicator = 0

	indicator = 1
	for i in range(N):
		if indicator == 1:
			if chi_list[-i] >= (max_value - delta_L):
				right_sigma = x_plot[-i]
				indicator = 0

	print("left_sigma",left_sigma,"(",abs(left_sigma - max_x),")")
	print("right_sigma",right_sigma,"(",abs(right_sigma - max_x),")")

	return left_sigma, right_sigma, max_x, max_value

#########################################################################

def plot_MLE(MLE,x_plot,max_value,max_x,right_sigma,left_sigma):
	plt.plot(x_plot,MLE)
	plt.axhline(y=max_value, color='red', linestyle='-', linewidth=1.4, label='BCG')
	plt.axhline(y=max_value-delta_L, color='red', linestyle='-', linewidth=1.4, label='BCG')

	plt.axvline(x=max_x, color='black', linestyle='-', linewidth=1.4, label='BCG')
	plt.axvline(x=right_sigma, color='black', linestyle='--', linewidth=1.4, label='BCG')
	plt.axvline(x=left_sigma, color='black', linestyle='--', linewidth=1.4, label='BCG')


	plt.show()



#########################################################################


MLE_A = mle(A,A_neg,A_poz,N)
MLE_B = mle(B,B_neg,B_poz,N)


MLE_sum = []

if str(operation) == "+":
	for i in range(N):
		MLE_sum = np.append(MLE_sum, MLE_A[i]+MLE_B[i])
elif str(operation) == "-":
	for i in range(N):
		MLE_sum = np.append(MLE_sum, (MLE_A[i]-MLE_B[i]))
elif str(operation) == "/":
	for i in range(N):
		MLE_sum = np.append(MLE_sum, -(MLE_A[i]/MLE_B[i]))
elif str(operation) == "*":
	for i in range(N):
		MLE_sum = np.append(MLE_sum, -(MLE_A[i]*MLE_B[i]))
else:
	print("\n#######################################")
	print("#        Invalid Operation!           #")
	print("#######################################\n")


x_plot_A = np.linspace(graphp_limits_A[0],graphp_limits_A[1],N)
x_plot_B = np.linspace(graphp_limits_B[0],graphp_limits_B[1],N)
x_plot_sum = np.linspace(graphp_limits_sum[0],graphp_limits_sum[1],N)

left_sigma_A, right_sigma_A, max_A, max_value_A = find_sigma(MLE_A,x_plot_A,N,delta_L)
left_sigma_B, right_sigma_B, max_B, max_value_B = find_sigma(MLE_B,x_plot_B,N,delta_L)
left_sigma_sum, right_sigma_sum, max_sum, max_value_sum = find_sigma(MLE_sum,x_plot_sum,N,delta_L)

#########################################################################

plot_MLE(MLE_sum,x_plot_sum,max_value_sum,max_sum,right_sigma_sum,left_sigma_sum)






