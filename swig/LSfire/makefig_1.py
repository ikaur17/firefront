# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 23:50:15 2013

@author: andrea
"""

import numpy as np
import matplotlib.pylab as plt
import os

try:
    LSpost_load
except:
    execfile('LSpost.py')
    print "...LSpost executed."
 
ft2m, m2ft = 0.3048, 1.0/0.3048
in2cm, cm2in = 2.54, 1.0/2.54
verbosity = 0

try: 
    ws_a = CreateWorkspace('LSout_s1_k0_U0670_I10_o0_b', verbose=verbosity)
    (it_a,t_a) = read_t(ws_a, verbose=verbosity)
except:
    (ws_a, it_a,t_a) = (None, None, None)
try: 
    ws_b = CreateWorkspace('LSout_s1_k0_U0670_I10_o1_b', verbose=verbosity)
    (it_b,t_b) = read_t(ws_b, verbose=verbosity)
except:
    (ws_b, it_b,t_b) = (None, None, None)
    
try: 
    ws_c = CreateWorkspace('LSout_s1_k1_U0670_I10_o0_b', verbose=verbosity)
    (it_c,t_c) = read_t(ws_c, verbose=verbosity)
except:
    (ws_c, it_c,t_c) = (None, None, None)
    
try: 
    ws_d = CreateWorkspace('LSout_s1_k1_U0670_I10_o1_b', verbose=verbosity)
    (it_d,t_d) = read_t(ws_d, verbose=verbosity)
except:
    (ws_d, it_d,t_d) = (None, None, None)

try: 
    ws_e = CreateWorkspace('LSout_s1_k3_U0670_I10_o0_b', verbose=verbosity)
    (it_e,t_e) = read_t(ws_e, verbose=verbosity)
except:
    (ws_e, it_e,t_e) = (None, None, None)
    
try: 
    ws_f = CreateWorkspace('LSout_s1_k3_U0670_I10_o1_b', verbose=verbosity)
    (it_f,t_f) = read_t(ws_f, verbose=verbosity)
except:
    (ws_f, it_f,t_f) = (None, None, None)
    
try:
    plt.close('all')
except:
    pass

print "time:", (t_a, t_b)
print "    :", (t_c, t_d)
print "    :", (t_e, t_f)


dim = (7,5) #(7*cm2in, 5*cm2in)

tt = """ k110: (368.496191, 228.11669, 228.11669)
         k112: (515.142634, 295.172969, 295.172969) """
L11 = [ 0, 50, 100, 150, 200, 225, 250, 300 ]


tt = """ k120: (127.010025, 83.559227, 83.350329)
         k122: (179.44344, 188.426057, 188.426057) """
L12 = [ 0, 25, 50, 75, 100, 125, 150]


tt = """ k210: (368.496191, 200.542145, 199.915451)
         k212: (440.566024, 387.297017, 364.736026) """
L21 = [ 0, 50, 100, 150, 198, 250, 300 ]


tt = """ k220: (127.010025, 66.638484, 66.638484)
         k222: (157.718041, 235.63702, 76.247795) """
L22 = [ 0, 25, 50,  65, 100, 150, 200 ]

subfig_a, subfig_c, subfig_e = '(a)', '(c)', '(e)'
subfig_b, subfig_d, subfig_f = '(b)', '(d)', '(f)'

x0, x1, y0, y1 = 0, 6000, 0, 6000 
fs = 14 #8
bbox = 'tight'

x_label = 'x [m]'
x_ticks = np.array([0,1,2,3,4,5])*1000*m2ft
x_ticks_label = ('0', '1000', '2000', '3000', '4000', '5000', '5000')

y_label = 'y [m]'
y_ticks = np.array([0,1,2,3,4,5])*1000*m2ft
y_ticks_label = ('0',  '1000',  '2000', '3000','4000','5000')

def common_makeup(fig, titolo, filename, cs=None, cloc=None):   
    ax = fig.gca()
    ax.axis([x0*m2ft, x1*m2ft, y0*m2ft, y1*m2ft]) ;
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_ticks_label, fontsize=fs)
    ax.set_xlabel(x_label, fontsize=fs)
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_ticks_label, fontsize=fs)    
    ax.set_ylabel(y_label, fontsize=fs)
    ax.set_title(titolo, fontsize=fs)
    if cloc is None:
        ax.clabel(cs, inline=1, fontsize=fs-1, fmt='%d')
    else:
        ax.clabel(cs, inline=1, fontsize=fs-1, fmt='%d', manual=cloc)
    plt.draw()
    fig.savefig(filename, bbox_inches=bbox)
    

###############################################################################
    
cpool = ( '#bd2309', '#bbb12d', '#1480fa', '#14fa2f', '#000000', #'#faf214', 
'#2edfea', '#ea2ec4', '#ea2e40', '#cdcdcd',
              '#577a4d', '#2e46c0', '#f59422', '#219774', '#8086d9' )
frames, show_frames, cpool = -1, True, None

Lstep20 = [20,40,60,80,100,120,140,160]#,180,200,220,240,260,280,300]
Lstep25 = [0,25,50,75,100,125,150,175,200,250] #,175,200,225,250,275,300]
try:
    La = [25, 50, 200, 100, 150]
    (ff, fig_a, cs1) = make_contour(ws_a, 'tburn_phieff', \
                         it=frames, show_frames=show_frames, L=Lstep25, dim=dim)
    common_makeup(fig_a, subfig_a, 'fig_a.pdf', cs1)
except:
    pass

try:
    Lb = [25, 50, 100, 150]
    (ff, fig_b, cs2) = make_contour(ws_b, 'tburn_phieff', \
               it=frames, show_frames=show_frames, L=Lb, dim=dim, cm=cpool)
    common_makeup(fig_b, subfig_b, 'fig_b.pdf', cs2)
except:
    pass

try:
    Lc = Lstep25
    (ff, fig_c, cs3) = make_contour(ws_c, 'tburn_phieff', \
               it=frames, show_frames=show_frames, L=Lc, dim=dim, cm=cpool)
    common_makeup(fig_c, subfig_c, 'fig_c.pdf', cs3)
except:
    pass
    
try:
    Ld = [25, 50, 100, 150]
    (ff, fig_d, cs4) = make_contour(ws_d, 'tburn_phieff', \
                         it=frames, show_frames=show_frames, L=Ld, dim=dim)
    common_makeup(fig_d, subfig_d, 'fig_d.pdf', cs4)
except:
    pass
    
try:
    Le = [0,10,20,25, 50,75,100, 200, 75, 100]
    (ff, fig_e, cs5) = make_contour(ws_e, 'tburn_phieff', \
               it=frames, show_frames=show_frames, L=Lstep25, dim=dim, cm=cpool)
    common_makeup(fig_e, subfig_e, 'fig_e.pdf', cs5)
except:
    pass
    
try:
    Lf = Lstep25
    (ff, fig_f, cs6) = make_contour(ws_f, 'tburn_phieff', \
                        it=frames, show_frames=show_frames, L=Lf, dim=dim)
    common_makeup(fig_f, subfig_f, 'fig_f.pdf', cs6)
except:
    pass

os.system('rm -f frame_tburn_phieff_001.png')
