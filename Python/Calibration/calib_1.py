# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 11:00:36 2018

@author: Nelis
"""
import scipy.io as spio
import numpy as np
import matplotlib.pyplot as plt
import math
import pylab as py

from readFile import read
from scipy import signal

counter = 0
res = 1001
N_runs = 1000
N_iter = 5
SNR = 10

shape = (res,1)
g_l = np.zeros(shape)
g_g = np.zeros(shape)
g_c = np.zeros(shape)
g_c_sampled = np.zeros(shape)

g_c_buff = np.zeros(shape)
g_c_new = np.zeros(shape)
g_buff = np.zeros(shape)

alpha = np.zeros(shape)
beta = np.zeros(shape)

shape = (91,1)
g_c_samples = np.zeros(shape)
g_l_samples = np.zeros(shape)
g_g_samples = np.zeros(shape)


src_ind = np.zeros(shape,dtype=np.int)
g_m = np.zeros(shape,dtype=np.int)
g_est = np.zeros(shape,dtype=np.int)
jitter = np.zeros(shape)

shape = (901,N_runs)
g_c_est = np.zeros(shape)

shape = (901,1)
g_c_est_std = np.zeros(shape)


l_range = np.linspace(-0.5,0.5,res)
l_range = l_range.reshape(len(l_range),1)
AF_1 = np.zeros((len(l_range),1))

shape = (901)
v = np.zeros(shape)
l_b = np.zeros(shape)

def import_beams():
    
    global beam_g1 
    beam_g1 = read('FF_11.ffe')
    beam_g2 = read('FF_inf_dense.ffe')
    
    py.plot(10**(beam_g1/10)/(max(10**(beam_g1/10))))
    py.plot(10**(beam_g2/10)/(max(10**(beam_g2/10))))
    
def gen_arr():
    
    global A_x
    
    N_ant = 20
    shape = (N_ant,1)
    A_x = np.ones(shape,dtype=complex)
    alpha = np.complex128(shape)
    beta = np.complex128(shape)
    
    lambda_0 = 0.3
    d_i = lambda_0/2
    k = 2*math.pi/lambda_0
    d_vec = np.arange(0,N_ant*d_i,d_i)
    
    l_0 = np.sin(20/180*np.pi)
    A_x = np.exp(-1j*k*d_vec*l_0)
    
#    for i in range(len(l_range)):
#        AF_1[i] = np.sum(np.sin(k*N_ant*d_vec*(1-l_range[i])/2)/(np.sin(k*d_vec*(1-l_range[i]))))

    for i in range(len(l_range)): 
        alpha = np.transpose(A_x)
        beta = np.exp(1j*k*d_vec*l_range[i]).reshape(N_ant,1)
        
        AF_1[i] =  np.dot(alpha,beta)
    
#    py.plot(l_range,AF_1)
    
def gen_beam():
    
    global l_b
    global counter
    global v
    
    v = np.linspace(0,1,901)
    l_b = 0.5*signal.sawtooth(2*np.pi*9*v)
    l_b[900] = 0.5
        
    for i in range(res):
        g_g[i] = 1 + -1*np.power(l_range[i], 2)
        
        if l_range[i] >= -0.45 and l_range[i] < -0.35:
            l_0 = -0.4
            g_l[i] = 1 - 200*np.power(l_range[i]-l_0, 2)
            counter = counter+1
                    
        if l_range[i] >= -0.35 and l_range[i] < -0.25:
            l_0 = -0.3
            g_l[i] = 1 - 200*np.power(l_range[i]-l_0, 2)

 
        if l_range[i] >= -0.25 and l_range[i] < -0.15:
            l_0 = -0.2
            g_l[i] = 1 - 200*np.power(l_range[i]-l_0, 2)
        

        if l_range[i] >= -0.15 and l_range[i] < -0.05:
            l_0 = -0.1
            g_l[i] = 1 - 200*np.power(l_range[i]-l_0, 2)    
        

        if l_range[i] >= -0.05 and l_range[i] < 0.05:
            l_0 = -0
            g_l[i] = 1 - 200*np.power(l_range[i]-l_0, 2)
        

        if l_range[i] >= 0.05 and l_range[i] < 0.15:
            l_0 = 0.1
            g_l[i] = 1 - 200*np.power(l_range[i]-l_0, 2)
        

        if l_range[i] >= 0.15 and l_range[i] < 0.25:
            l_0 = 0.2
            g_l[i] = 1 - 200*np.power(l_range[i]-l_0, 2)
        

        if l_range[i] >= 0.25 and l_range[i] < 0.35:
            l_0 = 0.3
            g_l[i] = 1 - 200*np.power(l_range[i]-l_0, 2) 
            
  
        if l_range[i] >= 0.35 and l_range[i] <= 0.45:
            l_0 = 0.4
            g_l[i] = 1 - 200*np.power(l_range[i]-l_0, 2)

        
        g_l[950] = 0.5
        g_c[i] = g_l[i] * g_g[i]
        
        
        
       
def MC_sim():
    
    ## define/state globals
    global g_m
    global jitter
    global g_est
    global g_l_est
    global g_g_est
    global g_c_est
    global alpha
    global beta
    global src_dstb

    ## Create regular pattern samples with regular source distribution
    src_dstb = np.linspace(50,950,91)

    ## Create random pattern samples with random source distribution
#    b_samples = 1*9
#    shape = (b_samples,1)
#    src_dstb = np.zeros(shape,dtype=np.int)
#    src_dstb = math.factorial(len(l_range))/math.factorial(len(l_range)-b_samples)

    for i in range(91):
        src_ind[i] = int(src_dstb[i])
    
    counter = 0  
    for i in range(50,951):
        if i in src_ind:
            g_c_sampled[i] = g_c[i]
            g_c_samples[counter] = g_c[i]
            g_l_samples[counter] = g_l[i]
            g_g_samples[counter] = g_g[i]
            
            counter += 1
            
    ## Add jitter to samples (forced data corruption)  
    for run_nr in range (N_runs):
        for i in range(91):
            jitter[i] = SNR + np.random.normal(1,1,1)
            
        g_m = g_c_samples * jitter
    
    ## Start multi-beam solution (iterative method)
        shape = (901,1)
        g_l_est = np.ones(shape)
        g_g_est = np.ones(shape)
    
        for iter in range(N_iter):
            # Remove global estimate from total estimation
            g_est = g_m.reshape(91)/g_g_est[src_ind-50].reshape(91)
            
            # Update calibration gain estimates for local model
            alpha = np.polyfit(l_b[src_ind-50].reshape(91), g_est, 2)
            g_l_est = np.polyval(alpha, l_b)
                        
            # Remove local estimate from total estimation
            g_est = g_m.reshape(91)/g_l_est[src_ind-50].reshape(91)
            
            # Update calibration gain estimates for global model
            beta = np.polyfit(l_range[src_ind].reshape(91), g_est, 2)
            g_g_est = np.polyval(beta, l_range[50:951])            
        
        g_c_est[:,run_nr] = g_l_est.reshape(901) * g_g_est.reshape(901)
                
def plot_stats():
    global g_c_est_std
    global buffer
    
    buffer = np.transpose(g_c_est/SNR)
    shape = (1,N_iter)
        
    for i in range(901):
        g_c_est_std[i] = np.std(buffer[:,i])
    
    
    image1 = plt.figure()
    py.plot(np.linspace(-0.5,0.5,901),g_c_est_std)
    plt.xlabel('l (direction cosine)')
    plt.ylabel('standard deviation')
    
    image1.savefig('im1.jpg', format = 'jpg',dpi=600)


def plot_beam():
    
#    py.plot(l_range,g_g) 
#    py.plot(l_range,g_l)

    image1 = plt.figure()
#
#    plt.plot(l_range[50:949],g_l[50:949],'--')
#    plt.plot(l_range[50:949],g_c[50:949])
#    plt.plot(l_range[50:949],g_g[50:949])
    

    plt.plot(l_range[50:949],g_c[50:949])
#    plt.plot(l_range,g_c_sampled,'rx')
    
    for i in range(50,949):
        if i in src_ind:
            plt.plot(l_range[i],g_c[i],'rx')
    
    plt.xlabel('l (direction cosine)')
    plt.ylabel('gain (normalised to unity)')

    image1.savefig('im1.jpg', format = 'jpg',dpi=600)


#import_beams()  
gen_arr()  
gen_beam()
MC_sim()
plot_stats()
#plot_beam()
