
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

#AF_full = np.ones(len(l_range),dtype=np.complex_)

sols = np.zeros(N_elements)

shape = (N_elements,1)
A_x = np.ones(shape,dtype=complex)
A_est = np.zeros(shape)
AF_samples = np.zeros(shape)

alpha = np.complex128(shape)
beta = np.complex128(shape)

shape = (len(l_range),1)
AF_full = np.complex128(shape)

shape = (N_elements,N_elements)
a_mat = np.zeros(shape)
a_inv = np.zeros(shape)

shape = (len(l_range),1)
AF_full = np.complex128(shape)


## Add jitter to A_vec
#for i in range(N_elements):
#   A_x[i] = 1 + np.random.normal(0,0.1,1)
   
## Construct constant A_vec   
#A_x[20] = 0
#A_x[10] = 0
#A_x[35] = 0   


## Calculate AF over l with jittered A_vec   
for i in range(len(l_range)): 
    alpha = np.transpose(A_x)
    beta = np.exp(1j*k*D_vec*l_range[i]).reshape(50,1)

    x = np.real(np.dot(alpha,beta)) + 1j*np.imag(np.dot(alpha,beta))
    AF_full[i] =  np.dot(alpha,beta)


### Take N samples from AF to solve A_x 
#for i in range(N_elements):
#    AF_samples[i] = AF_full[350+i]
#
### Set up A_matrix
#for i in range(N_elements):
#    for j in range(N_elements):
#        a_mat[i][j] = np.exp(1j*k*j*d_i*l_range[350+i])
#        
## solve x_vec
#a_inv = np.linalg.inv(a_mat)
#x_vec = np.matmul(a_inv,AF_samples)

    AF_full[i] =  np.dot(alpha,beta)



## Evaluate Cramer
#for i in range(N_elements):
#    sols[i] = AF_full[i]


    
#py.plot(AF_samples)
#py.plot(20*np.log10(AF_full))