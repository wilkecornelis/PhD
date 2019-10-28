# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 11:01:43 2018

@author: Nelis
"""
import numpy as np

def read(filename):
            
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    data = lines[15:len(lines)-1]
    rows = [line.split() for line in data]
    
    extract = np.zeros((len(rows),1))    
    
    for i in range(len(rows)):
        buffer = rows[i]
        extract[i] = buffer[8]
    
    return extract
    