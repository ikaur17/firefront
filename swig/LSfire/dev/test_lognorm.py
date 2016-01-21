# -*- coding: utf-8 -*-
"""
Created on Sat Jul  6 12:06:59 2013

@author: andrea
"""

Nquad = 100
nfquad_x = '../quadlib/quad_hermite_'+str(Nquad)+'_x.txt'
nfquad_w = '../quadlib/quad_hermite_'+str(Nquad)+'_w.txt'

x_q, w_q = [], []

f_x = open(nfquad_x, 'r')
for x in f_x.readlines():
    x_q.append(float(x))
f_x.close()

f_w = open(nfquad_w, 'r')
for w in f_w.readlines():
    w_q.append(float(w))
f_w.close()

def g0(csi, eta, sigma, mu, s):
    f = 2 * pi**1.5 * sigma**2
    return exp(-(csi**2+eta**2)/(2*sigma**2)-mu**2/(2*s**2)) / f
    
def g1(csi, sigma, mu, s):
    
    def f_q(w):
        t1 = -(exp(2*sqrt(2)*s*w)-2*csi*exp(sqrt(2)*s*w))/(2*sigma**2)
        t2 = sqrt(2)*mu*w/s
        return exp(t1+t2)

    S = 0
    for i in range(Nquad):
        S += w_q[i] * f_q(x_q[i])
    return S
    
def g(csi, eta, sigma, mu, s):
    G0, G1 = g0(csi, eta, sigma, mu, s), g1(csi, sigma, mu, s)
    return (G0*G1, G0, G1)