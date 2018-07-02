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

res = 1001
N_runs = 1000
N_iter = 5
SNR = 10

shape = (res,1)
g_l = np.zeros(shape)
g_g = np.zeros(shape)
g_c = np.zeros(shape)

g_c_buff = np.zeros(shape)
g_c_new = np.zeros(shape)
g_buff = np.zeros(shape)

alpha = np.zeros(shape)
beta = np.zeros(shape)

shape = (91,1)
src_dstb = np.zeros(shape,dtype=np.int)
g_c_samples = np.zeros(shape)
g_l_samples = np.zeros(shape)
g_g_samples = np.zeros(shape)


src_ind = np.zeros(shape,dtype=np.int)
g_m = np.zeros(shape,dtype=np.int)
g_est = np.zeros(shape,dtype=np.int)
jitter = np.zeros(shape)

shape = (1001,N_runs)
g_c_est = np.zeros(shape)

shape = (1001,1)
g_c_est_std = np.zeros(shape)


l_range = np.linspace(-0.5,0.5,res)

def gen_beam():
    
    for i in range(res):
        g_g[i] = 1 + -1*np.power(l_range[i], 2)
        
        if l_range[i] >= -0.45 and l_range[i] < -0.35:
            l_0 = -0.4
            g_l[i] = 1 - 200*np.power(l_range[i]-l_0, 2)
            
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

    ## Create pattern samples with source distribution
    src_dstb = np.linspace(1,len(g_c)-1,91)

    for i in range(91):
        src_ind[i] = int(src_dstb[i])
    
    counter = 0  
    for i in range(1001):
        if i in src_ind:
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
        shape = (1001,1)
        g_l_est = np.ones(shape)
        g_g_est = np.ones(shape)
    
        for iter in range(N_iter):
            # Remove global estimate from total estimation
            g_est = g_m.reshape(91)/g_g_est[src_ind].reshape(91)
            
            # Update calibration gain estimates for local model
            alpha = np.polyfit(g_l_samples.reshape(91), g_est, 2)
            g_l_est = np.polyval(alpha, g_l)
                        
            # Remove local estimate from total estimation
            g_est = g_m.reshape(91)/g_l_est[src_ind].reshape(91)
            
            # Update calibration gain estimates for global model
            beta = np.polyfit(l_range[src_ind].reshape(91), g_est, 2)
            g_g_est = np.polyval(beta, l_range)            
        
        g_c_est[:,run_nr] = g_l_est.reshape(1001) * g_g_est.reshape(1001)
        
        
def plot_stats():
    global g_c_est_std
    global buffer
    
    buffer = np.transpose(g_c_est/SNR)
    
    for i in range(1001):
        g_c_est_std[i] = np.std(buffer[:,i])
    
    
    py.plot(g_c_est_std)
    py.plot(g_c/50)
    
def plot_beam():
    
#    py.plot(l_range,g_g) 
#    py.plot(l_range,g_l)
    py.plot(l_range,g_c)


    
gen_beam()
MC_sim()
plot_stats()
#plot_beam()
