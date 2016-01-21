# -*- coding: utf-8 -*-
"""
Created on Fri May  3 21:02:26 2013
@author: andrea
"""

from pylab import *
from scipy.sparse import spdiags


h0 = 30 ; D = (h0**2)/4.0 ;
dr = 10 ; Nr = 100 ;
dt = 1 ; Nit = 20 ;

dt = 1 ; dr = 10; Nr = 1414 ; Nit = 10 ; #10, 1414, 15

flag_method = 21
#flag_adapt = 1
GRAD_P_MIN = 1e-12

tmax = dt * Nit

def coeff():

    global r, Nr, p0, a, b, c, A
    global Dr, ri, DR

    (r0, r1) = (0.0, (Nr-1)*dr)
    
    r = arange(r0, r1+dr, dr)
    Nr = size(r)
    
    #print "........Nr: "+str(Nr)+", Nit: "+str(Nit)+" r1): "+str(r1)
    
    Dr = r[1:] - r[:-1]
    ri = (r[:-1] + r[1:]) / 2.0
    DR = r * ( hstack([ri, ri[-1]+Dr[-1]]) - hstack([ri[0]-Dr[0], ri]) )
    
    sx = hstack([ nan, ri/Dr ]) / DR
    dx = hstack([ nan, ri[1:]/Dr[1:], (ri[-1] + Dr[-1])/Dr[-1] ]) / DR
    cc = -sx-dx ; 
    cc[0] = -2.0 / Dr[0]**2
    dx[0] = 2.0 / Dr[0]**2
    
    (c, a, b) = (sp.hstack([sx[1:], nan]), cc, hstack([nan, dx[0:-1]]) )
    
    if flag_method == 10 or flag_method == 11:
        A = spdiags(vstack([c, a, b]), [-1, 0, 1], Nr, Nr)
    else:
        A = []

def disp():
    if not ((100.*it)/Nit)%10:
        print "it #" + str(it) + ", t=" + str(t) + " {Nr: " + str(Nr) + "} > max(p)=" + str(max(p))
        
def f(t, p):
    return D * A * p
   
###############################################################################
print ""

coeff()
#print "r[i] = ", r, "\n\nDr[i] = ", Dr, "\n\nri[i] = ", ri, "\n\nDR[i] = ", DR
#print "\n\nc[i] = ", c[:-1], "\n\na[i] = ", a, "\n\nb[i] = ", b[1:]

p0 = zeros(Nr); p0[0] = 1.0 / (pi*(r[1]/2.0)**2) / 2.0

#pzero=p0

t0 = 0 ; it = 0 ; 
(t, t_old, p_old) = (0.0, 0.0, p0)
NR_MIN, NR_MAX = 10, 1000000

# flag_method = 10 : Runge Kutta 4/5
# flag_method = 11 : Runge Kutta 4/5 + corrector
# flag_method = 20 : implicito 
# flag_method = 21 : implicito (ottimizzato)

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
        p = solve(eye(Nr)-dt*D*A, p_old )
        
    elif flag_method == 21:
        p = thomas([-dt*D*c, ones(Nr)-dt*D*a, -dt*D*b], p_old)
        
    if flag_adapt > 0 and it<Nit:
        grad_p = abs((p[1:]-p[:-1]) / (r[1:]-r[:-1])) 
        last_grad = grad_p[-1]
        Nr_old = Nr
        if flag_adapt == 2:
            Nr = len(grad_p[grad_p>GRAD_P_MIN]) + 1
            Nr = min( max(Nr, NR_MIN), Nr_old)
            #print "..(2) Nr_old, Nr:", Nr_old, Nr
        if flag_adapt >= 1 and last_grad > GRAD_P_MIN :
            Nr = min( max(int(ceil(Nr*1.05)), NR_MIN), NR_MAX)
            #print "..(1) Nr_old, Nr:", Nr_old, Nr
        if Nr != Nr_old:
            r1 = (Nr-1)*dr
            coeff()
            if Nr > Nr_old:
                p_old = hstack([p_old, zeros(Nr-Nr_old)])
            else:
                p_old = p_old[:Nr]
            disp()
            continue
    
    t_old, p_old = t, p
    t += dt
    
    disp()
    
    if isnan(max(p)): break

###############################################################################

def p_exact(r, t):
    return exp(-r**2/(4*D*t)) / (4*pi*D*t)
    
close('all')
pexact = p_exact(r, t)
plot(r, 2*pi*r*pexact, 'r', r, 2*pi*r*p, '.-b')
grid('on')
print "...max_e:", max(abs(p-pexact)/p*100), "%,", max(abs(r*p-r*pexact)/(r*p)*100), "%"