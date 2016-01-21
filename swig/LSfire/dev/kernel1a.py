# -*- coding: utf-8 -*-
"""
Created on Fri May  3 21:02:26 2013

@author: andrea
"""

from pylab import *

import scipy as sp
from scipy.sparse import spdiags

sigma = .001
(r0, r1) = (0.00000001, 10.0)
N = 100

tmax = 1000 ; dt = 1

p0 = sp.zeros(N)
p0[0] = 1

r = linspace(r0, r1, N)
#r = (logspace(r0, r1, N)-10**r0) / (10**r1-10**r0)*r1
dr = (r[1:]-r[:-1])
ri = (r[:-1]+r[1:])/2.0

flag_A = 2
if flag_A == 1:
    cc = -2*ones(N) 
    dx = hstack([ri/r[:-1], 0]) ; dx[1] = -cc[0] #2*dx[1]
    sx = hstack([ri/r[1:], 0.0])
    cc = cc / dr[0]**2
    dx = dx / dr[0]**2
    sx = sx / dr[0]**2
    A = spdiags(sp.vstack([sx, cc, dx]), [-1, 0, 1], N, N)
else:
    Dr = r * ( hstack([ri,ri[-1]+dr[-1]]) - hstack([ri[0]-dr[0],ri]) )
    sx = sp.hstack([ (ri[0]-dr[0])/dr[0], (ri/dr) ]) / Dr
    dx = sp.hstack([ (ri/dr), (ri[-1]+dr[-1])/dr[-1] ]) / Dr
    cc = -sx-dx
    cc[0]=-2.0/dr[0]**2 ; dx[0] = 2.0/dr[0]**2
    #cfr = ( (-2.0/dr[0]**2, 2.0/dr[0]**2), (cc[0], dx[0]+sx[0]) )
    #dx[0] += sx[0]
    A = spdiags(sp.vstack([sp.hstack([sx[1:],nan]), cc, sp.hstack([nan,dx[0:-1]])]), [-1, 0, 1], N, N)




def f(t, p):
    c = sigma#/ (2*t+1e-6) 
    return c * A * p

t0 = 0; it = 0 ; 
t, t_old, p_old = 0.0, 0.0, p0
   
while t < tmax:
    
    it += 1
        
    k1 = f(t_old, p_old)
    k2 = f(t_old+dt/2.0, p_old+k1/2.0)
    k3 = f(t_old+dt/2.0, p_old+k2/2.0)
    k4 = f(t_old+dt, p_old+k3)
    
    t += dt
    p = p_old + (k1 + 2*k2 + 2*k3 + k4)*dt/6.0
    
    t_old, p_old = t, p
    
    print "it #" + str(it) + ", t=" + str(t) + " > max(p)=" +\
        str(max(p))
        
    if isnan(max(p)): break
close('all')
#plot(r, p0,'r', r, p, '.-b')
plot(r, p, '.-b')
    #fig = plot(p)
    #show(fig)
    #pause(1)