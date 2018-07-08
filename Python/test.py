# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 13:24:12 2018

@author: Nelis
"""


import numpy as np
import math
import pylab as py

lambda_0 = 0.3
d_i = lambda_0/2
N_elements = 50
k = 2*math.pi/lambda_0

l=0.5
theta_range = np.arange(-math.pi/2,math.pi/2,0.005)
l_range  = np.sin(theta_range)
A_vec = np.ones(N_elements)
D_vec = np.arange(0,N_elements*d_i,d_i)
AF_full = np.ones(len(l_range))
sols = np.zeros(N_elements)

shape = (N_elements,1)
A_x = np.zeros(shape)
A_est = np.zeros(shape)

## Add jitter to A_vec
for i in range(N_elements):
   A_x[i] = 1 + np.random.normal(0,0.1,1)

## Calculate AF over l with jittered A_vec   

for i in range(len(l_range)): 
    alpha = np.transpose(A_x)
    beta = np.exp(1j*k*D_vec*l_range[i]).reshape(50,1)
    AF_full[i] =  abs(np.dot(alpha,beta))



## Evaluate Cramer
for i in range(N_elements):
    sols[i] = AF_full[i]

    



 