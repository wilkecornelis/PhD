# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 11:40:43 2018

@author: Nelis
"""

import numpy as np
import random
import matplotlib.pyplot as plt
import math
import pylab

N_sources = 5
N_ant = 7

x = 512
y = 512

a = [[0]* y for _ in range(x)]

def column(matrix, i):
    return [row[i] for row in matrix]
            
def gen_array_coords():

    
    y_m = np.array([25.1,90.28,3.99,-21.6,-38.27,-61.6,-87.99])
    x_m = np.array([-9.1,26.38,26.89,25.49,-2.59,-79.7,75.75])
    
    global x_lambda
    global y_lambda

    x_lambda = [0]* N_ant
    y_lambda = [0]* N_ant                      
    
    for i in range(N_ant):
        x_lambda[i] = x_m[i]/0.2
        y_lambda[i] = y_m[i]/0.2        
                
                
    plt.scatter(x_lambda,y_lambda)
    plt.xlabel('X_Wavelengths')
    plt.ylabel('Y_Wavelengths')
    plt.show()           
    
def gen_uv_coords():
    
    H = 0
    decl = 90*math.pi/180     
    
    uv_mat_len = math.factorial(N_ant)/math.factorial(N_ant - 2) ## Calculate number of U and V coords using permutations
                                                            
    u = [0]* int(uv_mat_len)
    v = [0]* int(uv_mat_len)   
    c = 0
    
    for i in range(N_ant):
        for j in range(N_ant):
            if i != j:
                Dx = x_lambda[i]-x_lambda[j]
                Dy = y_lambda[i]-y_lambda[j]
                
                u[c] = math.sin(H)*Dx + math.cos(H)*Dy
                v[c] = -1*math.sin(decl)*math.cos(H)*Dx + math.sin(decl)*math.cos(H)*Dy 
                
                c = c+1
                
    plt.scatter(u,v)
    plt.show()

def gen_dyn_uv_coords():
    
    decl = 50*math.pi/180       ## declination of source in radions
    h_angle = 6                 ## observation range in hours
    h_interval = 15             ## sample interval in minutes
    h_range = (h_angle/6)*90    ## calculate observation range in degrees for given observation duration
    
    N_intervals = h_range/((h_interval/60)*15)      ## number of UV sample intervals for given observation in hour angles
    interval_angle = (h_interval/60)*15
                                                                                
    uv_mat_len = math.factorial(N_ant)/math.factorial(N_ant - 2) ## Calculate number of U and V coords using permutations nPr of baselines
    
    global u_d
    global v_d

    global u_f
    global v_f    
                                                            
    u_d = [0]* int(uv_mat_len*int(N_intervals))
    v_d = [0]* int(uv_mat_len*int(N_intervals))   
    
    u_f = [0]* int(uv_mat_len*int(N_intervals))
    v_f = [0]* int(uv_mat_len*int(N_intervals))  
    
    global c1
    c1 = 0
    
    for k in range(int(N_intervals)): 
        H = k*interval_angle*math.pi/180
        
        for i in range(N_ant):
            for j in range(N_ant):
                if i != j:
                    Dx = x_lambda[i]-x_lambda[j]
                    Dy = y_lambda[i]-y_lambda[j]
                    
                    u_d[c1] = math.sin(H)*Dx + math.cos(H)*Dy
                    v_d[c1] = -1*math.sin(decl)*math.cos(H)*Dx + math.sin(decl)*math.cos(H)*Dy 
                       
                    u_f[c1] = 300e6/(math.sin(H)*Dx + math.cos(H)*Dy)
                    v_f[c1] = 300e6/(-1*math.sin(decl)*math.cos(H)*Dx + math.sin(decl)*math.cos(H)*Dy)
                    
                    c1 = c1+1 
            
    plt.scatter(u_d,v_d)
    plt.xlabel('v')
    plt.ylabel('u')
    plt.show()
#    plt.scatter(u_d,v_d)
#    plt.show()    

def gen_dirty_beam():
    
    print(len(u_f))
    
    nrows = 16000
    ncols = 16000
    
    arr = np.zeros((nrows, ncols))
    
    for i in range (c1):
        if v_d[i] > 0:
            if u_d[i] > 0:
                arr[int(v_f[i])][int(u_f[i])] = 1 
    
    B = np.fft.ifft2(arr)
    
#    pylab.imshow(np.abs(B))
#    pylab.show()
    
    l = np.arange(nrows)
    m = np.arange(ncols)
    L, M = np.meshgrid(l, m)
    plt.contourf(L, M, np.abs(B))
    plt.show()        

def cal_visibility():
    A = np.exp(-2*np.pi*1j*())
    
                 
def gen_sky():
    
    source_x = random.randint(20,480)
    source_y = random.randint(20,480)
    
    
    print(source_x)   
    print(source_y)              
    
    global a

    for i in range(y):
        for j in range(x):
            a[i][j] = 0
         

    for i in range(N_sources):
        source_x = random.randint(20,480)
        source_y = random.randint(20,480) 
            
        a[source_x][source_y] = 1
    
def skyplot():
    n = 512
    l = np.arange(n)
    m = np.arange(n)
    L, M = np.meshgrid(l, m)
    plt.contourf(L, M, a)
    plt.show()
    
gen_array_coords()
gen_dyn_uv_coords()
#gen_dirty_beam()

#cal_visibility()


#gen_uv_coords()

#gen_uv_coverage()    
#gensky()    
#skyplot()




