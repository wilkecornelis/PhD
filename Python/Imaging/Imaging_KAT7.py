# -*- coding: utf-8 -*-
"""
Created on Wed May 16 14:08:21 2018

@author: Nelis
"""

import scipy.io as spio
import numpy as np
import matplotlib.pyplot as plt
import math

from matplotlib import ticker, cm

"""
Create variables and initialise
"""

mat = spio.loadmat('antPos.mat')
ant_pos = mat['posLBA']
f_s = 1000e6
c = 300e6
lambda_s = c/f_s
N_ant = 7

theta_s = 0*np.pi/180 
phi_s = 0*np.pi/180 

B_s = [[0]* 500 for _ in range(500)]
l_range = np.arange(-0.1,0.1,0.004)
m_range = np.arange(-0.1,0.1,0.004)

shape = (N_ant,1)
b = np.zeros(shape, np.complex128) 
a1 = np.zeros(shape, np.complex128)
a2 = np.zeros(shape, np.complex128)

shape = (1,N_ant)
b_T = np.zeros(shape)
a_T = np.zeros(shape, np.complex128)

shape = (500,500)
sky1 = np.zeros(shape)
sky2 = np.zeros(shape)

l_s = math.sin(theta_s)*math.cos(phi_s)
m_s = math.sin(theta_s)*math.sin(phi_s)

shape  = (N_ant,N_ant)
R_s1 = np.zeros(shape, np.complex128)
R_s2 = np.zeros(shape, np.complex128)

"""
Functions and methods
"""
    
def cal_covariance_mat():
   
    global ant_pos
    global R_s1    
    
    global y_m
    global x_m

    y_m = np.array([25.1,90.28,3.99,-21.6,-38.27,-61.6,-87.99])
    x_m = np.array([-9.1,26.38,26.89,25.49,-2.59,-79.7,75.75])
    

#    xl = l_s*ant_pos[:,0]
#    ym = m_s*ant_pos[:,1]
#    a = np.exp(-2*np.pi*1j*f_s*(xl+ym)/c)
    
    for i in range(N_ant):
        xl = l_s*x_m[i]
        ym = m_s*y_m[i]    
        g1 = 1
#        g1 = np.random.normal(1,0.01,1)
        a1[i] = g1*np.exp(-2*np.pi*1j*f_s*(xl+ym)/c) 

    a_T1 = np.conj(a1).T
    R_s1 = a1*a_T1

    
def cal_sky_mat():
    
    global sky1

    for l in range(500):
        for m in range(500):
                        
            xl = x_m[:]*l_range[l]
            ym = y_m[:]*m_range[m]         
            b = np.exp(-2*np.pi*1j*f_s*(xl+ym)/c)        
    
            b_T = np.conj(b).T
            B_s11 = np.matmul(b_T,R_s1)                 
            B_s21 = np.matmul(B_s11,b)
            
            sky1[l][m] = abs(B_s21)
     
def plot_sky():
           
    L, M = np.meshgrid(l_range, m_range)   
    
    plot1 = 10*np.log10(sky1)

    image1 = plt.figure()
        
    plt.contourf(L, M, sky1, cmap=cm.hsv)
    plt.xlabel('l')
    plt.ylabel('m')
    
    plt.colorbar()
    
    plt.show()
    image1.savefig('im1.jpg', format = 'jpg',dpi=600)


"""
Call functions and methods
"""

cal_covariance_mat()
cal_sky_mat()
plot_sky()

              