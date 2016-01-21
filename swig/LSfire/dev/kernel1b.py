# -*- coding: utf-8 -*-
"""
Created on Fri May  3 21:02:26 2013

@author: andrea
"""

from pylab import *

#import scipy as sp
from scipy.sparse import spdiags

D = 0.001 ; #D = 225.0 
dt = 1 ; tmax = dt * 300
(r0, r1) = (0.01, 9*(4*D)*tmax)#(9*h0**2)*tmax / 10)
N = 100

p0 = zeros(N); p0[0] = 1

r = linspace(r0, r1, N)
dr = (r[1:]-r[:-1])
ri = (r[:-1]+r[1:])/2.0

Dr = r * ( hstack([ri,ri[-1]+dr[-1]]) - hstack([ri[0]-dr[0],ri]) )
sx = hstack([ (ri[0]-dr[0])/dr[0], (ri/dr) ]) / Dr
dx = hstack([ (ri/dr), (ri[-1]+dr[-1])/dr[-1] ]) / Dr
cc = -sx-dx 
cc[0]=-2.0/dr[0]**2 ; dx[0] = 2.0/dr[0]**2
A = spdiags(vstack([sp.hstack([sx[1:],nan]), cc, hstack([nan,dx[0:-1]])]), [-1, 0, 1], N, N)

def exact_solution(r, t):
    return exp(-r**2/(4*D*t)) / (4*pi*D*t)

def f(t, p):
    return D * A * p

t0 = 0; it = 0 ; 
t, t_old, p_old = 0.0, 0.0, p0
   
while t < 100:#tmax:
    
    it += 1
        
    k1 = f(t_old, p_old)
    k2 = f(t_old+dt/2.0, p_old+k1/2.0)
    k3 = f(t_old+dt/2.0, p_old+k2/2.0)
    k4 = f(t_old+dt, p_old+k3)
    
    t += dt
    p = p_old + (k1 + 2*k2 + 2*k3 + k4)*dt/6.0
    
    t_old, p_old = t, p
    
    print "it #" + str(it) + ", t=" + str(t) + " > max(p)=" + str(max(p))
    
    if isnan(max(p)): break

close('all')
pexact = exact_solution(r, t)
plot(r, pexact, 'r', r, p, '.-b')
print "(sol_exact/sol_approx)[0] = ", pexact[0]/p[0]
    #fig = plot(p)
    #show(fig)
    #pause(1)