# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 22:06:43 2013

@author: andrea
"""

# https://github.com/johannesgerer/jburkardt-f/tree/master/gen_laguerre_rule

def fun(x, csi, t, h=2, D=1, lam=1):
    return exp(-(csi-x)**2/(2.0*D*t)-(x/lam)**h+x)
    
def get_gl(h, N):
    
    A, B = 0, 1    
    nomefile = 'gl_alpha' + str(h-1) + '_A' + str(A) + '_B' + str(B) + '_N' + str(N)
    
    x, w = [], []
    
    file_x = open(nomefile+'_x.txt','r')
    for line in file_x.readlines():
        x.append(float(line))
    file_x.close()
    
    file_w = open(nomefile+'_w.txt','r')
    for line in file_w.readlines():
        w.append(float(line))
    file_w.close() 
    
    return (x, w)
    

def get_integral(csi, t, N, h=2, D=1, lam=1):
    (x, w) = get_gl(h, N)
    I = 0
    for i in range(N):
        I += fun(x[i], csi, t, h, D, lam) * w[i]
        
    return I
    
# get_integral(.001, 1000, 100, 2, 1, 1000)