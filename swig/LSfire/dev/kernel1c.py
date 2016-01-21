# -*- coding: utf-8 -*-
"""
Created on Fri May  3 21:02:26 2013
@author: andrea
"""

from pylab import *
from scipy.sparse import spdiags

D = 225 # 0.001 ; #D = 225.0 
dr = 1 ; dt = 10 ; Nit = 100 ; # 0.1, 0.01, 1000
(r0, r1) = (0, 1000.) # 1e-6, 0.8

tmax = dt * Nit
Nz = 100
flag_method = 21

def coeff():
    global r, N, p0, a, b, c, A
    r = arange(r0, r1+dr, dr)
    r = hstack([ linspace(r[0], r[2], 3), r[3:] ]) 
    N = size(r)
    
    p0 = zeros(N); p0[0] = 1.0 / (pi*(r[1]/2.0)**2) / 2.0
    
    print "...N: "+str(N)+", Nit: "+str(Nit)+" (r0, r1): ("+str(r0)+", "+str(r1)+")"
    
    Dr = (r[1:]-r[:-1])
    ri = (r[:-1]+r[1:])/2.0
    DR = r * ( hstack([ri,ri[-1]+Dr[-1]]) - hstack([ri[0]-Dr[0],ri]) )
    
    sx = hstack([ nan, (ri/Dr) ]) / DR
    dx = hstack([ nan, (ri[1:]/Dr[1:]), (ri[-1]+Dr[-1])/Dr[-1] ]) / DR
    cc = -sx-dx ; cc[0]=-2.0/Dr[0]**2 ; dx[0] = 2.0/Dr[0]**2
    
    (c, a, b) = (sp.hstack([sx[1:],nan]), cc, hstack([nan,dx[0:-1]]) )
    
    if flag_method == 10 or flag_method == 11:
        A = spdiags(vstack([c, a, b]), [-1, 0, 1], N, N)
    else:
        A = []

def f(t, p):
    return D * A * p

print ""
coeff()
t0 = 0; it = 0 ; 
t, t_old, p_old = 0.0, 0.0, p0
   
for it in range(1, Nit+1):
    
    if flag_method == 10 or flag_method == 11:
        k1 = f(t_old, p_old)
        k2 = f(t_old+dt/2.0, p_old+k1/2.0)
        k3 = f(t_old+dt/2.0, p_old+k2/2.0)
        k4 = f(t_old+dt, p_old+k3)
        p = p_old + (k1 + 2*k2 + 2*k3 + k4)*dt/6.0
        if flag_method == 11:
            p = (p_old + p)/2.0 + f(t, p)*dt/2.0 # corrector
            
    elif flag_method == 20:
        p = solve(eye(N)-dt*D*A, p_old )
        
    elif flag_method == 21:
        p = thomas([-dt*D*c, ones(N)-dt*D*a, -dt*D*b], p_old)
        
    if abs((p[-1]-p[-2])/(r[-1]-r[-2])) > 1e-12:
        r1 = r1*1.05
        Nold = N
        coeff()
        p_old = hstack([p_old, zeros(N-Nold)]) 
        continue
    
    t_old, p_old = t, p
    t += dt
    
    if not ((100.*it)/Nit)%10: 
        print "it #" + str(it) + ", t=" + str(t) + " > max(p)=" + str(max(p))
    
    if isnan(max(p)): break

def p_exact(r, t):
    return exp(-r**2/(4*D*t)) / (4*pi*D*t)
    
close('all')
pexact = p_exact(r, t)
plot(r, 2*pi*r*pexact, 'r', r, 2*pi*r*p, '.-b')
grid('on')
print "...max_e:", max(abs(p-pexact)/p*100), "%,", max(abs(r*p-r*pexact)/(r*p)*100), "%"
    #fig = plot(p)
    #show(fig)\n
    #pause(1)