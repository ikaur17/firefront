# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 08:23:24 2013

@author: andrea

usage: import LSpost
       ws = LSpost.CreateWorkspace('LS')
       LSpost.make_movie_contour(ws, 'phi', range(1,8), arange(0,1,0.2))
       
or: runfile('LSpost.py')
    ws = CreateWorkspace('LSout01')
    ros = read_array(ws, 'ros', 10)
    make_movie_contour(ws, 'phi', range(1,11), arange(0,1,0.2) )
"""

import numpy as np
import matplotlib.pylab as plt
import platform
import os

from numpy import inf

LSpost_load = True

if platform.node() == 'undy-VII':
    sim_dir  = '/media/CIRAM/LSfire_sim/'
else:
    sim_dir  = '/home/andrea/CIRAM/LSfire_sim/'
    
L1 = np.arange(0, 1.1, 0.1)
L2 = np.arange(0, 1.2, 0.2)
L5 = [0.1, 0.5, 1.0]

col_grigio30 = '#B2B2B2'

class Bunch(object):
    """ If you have a dictionary d and want to read its values 
         with the syntax x.foo instead of d['foo'], just do:
             x = Bunch(d)
        [http://stackoverflow.com/questions/2597278/python-load-variables-in-a-dict-into-namespace]
    """
    def __init__(self, adict):
        self.__dict__.update(adict)
        
    
###############################################################################
    
       
def CreateWorkspace(filename_ws, path_ws=None, verbose=True):
    
    filename_data = 'data'
    
    filename_data_000 = filename_ws + '_' + filename_data + '_000.txt'
    filename_data_end = filename_ws + '_' + filename_data + '_end.txt'
    
    if path_ws:
        if path_ws[-1] != '/': path_ws += '/'
        try:
            f = open(path_ws+filename_data_000, 'r')
            f.close()
        except:
            print "No way, bebe."
            return None
    else:
        try_path_ws1 = './'
        try_path_ws2 = './' + filename_ws + '/'
        try:
            f = open(try_path_ws1+filename_data_000, 'r')
            path_ws = try_path_ws1
            f.close()
        except:
            try:
                f = open(try_path_ws2+filename_data_000, 'r')
                path_ws = try_path_ws2
                f.close()
            except:
                print "No way, bebe."
                return None                        
                
    class InjectInNamespace():

        """ Inject variables from filename """
        
        from numpy import inf
        
        def __init__(self):
            if verbose:
                print "...namespace for '{}' injected from path '{}'.".format(\
                filename_ws, path_ws)
                
        f = open(path_ws+filename_data_000, 'r')
        for line in f.readlines(): exec(line)
        f.close()
        
        f = open(path_ws+filename_data_end, 'r')
        for line in f.readlines(): exec(line)
        f.close()

    ws = InjectInNamespace()
    
    ws.path_ws = path_ws
    
    # iterations saved
    if ws.save_each_iteration > 0:
        step = ws.save_each_iteration
        ws.Niters = range(step, ws.Niter, step)
        if ws.Niters[0] != 1:
            ws.Niters = [1] + ws.Niters
        if ws.Niters[-1] < ws.Niter: 
            ws.Niters = ws.Niters + [ws.Niter]
    else:
        ws.Niters = [ws.Niter]
    
    (ws.xx, ws.yy, ws.n, ws.m, \
     ws.xxn, ws.yyn, ws.nn, ws.mn, \
     ws.xxc, ws.yyc, ws.nc, ws.mc) = make_grid(ws, 0)
     
    if verbose:
        print "...Niter is: {}.".format(ws.Niter)

    return ws

    
###############################################################################


def make_grid(ws, verbose=1):
    
    # NB: I valori calcolati sono riferiti ai centroidi della griglia
    #    (non ai nodi)
    
    # Queste sono le coordinate dei bordi delle celle
    if ws.flag_saveformat_ghostboxes:
        xxn = np.linspace(ws.x_0_g, ws.x_1_g, ws.Nx_g+1) ;
        yyn = np.linspace(ws.y_0_g, ws.y_1_g, ws.Ny_g+1) ;
    else:
        xxn = np.linspace(ws.x_0, ws.x_1, ws.Nx+1) ;
        yyn = np.linspace(ws.y_0, ws.y_1, ws.Ny+1) ;
        
    # Queste sono le coordinate dei centroidi
    xx, yy = xxn[0:-1]+ws.dx/2.0, yyn[0:-1]+ws.dy/2.0
    
    # Queste sono le coordinate della griglia croppata (senza ghost cells)
    i0, i1, j0, j1 = ws.i0_f, ws.i1_f, ws.j0_f, ws.j1_f
    xxc, yyc = xx[i0:i1+1], yy[j0:j1+1]
    
    n, m = np.size(xx), np.size(yy)
    nc, mc = np.size(xxc), np.size(yyc)
    nn, mn = np.size(xxn), np.size(yyn)
        
    if verbose:
        print "\n...generated ordinary grid:  "+str(n)+"x"+str(m)+"."
        print "\n...generated cropped grid:   "+str(nc)+"x"+str(mc)+"."
        print "\n...generated interface grid: "+str(nn)+"x"+str(mn)+".\n"
    
    return (xx, yy, n, m, xxn, yyn, nn, mn, xxc, yyc, nc, mc)
    

###############################################################################
  
def DisplayInfo(ws):
    
    print ""
    print " Grid basic info for workspace '"+ws.filename_ws+"':"
    print "   (x_0, x_1) = (", ws.x_0, ",", ws.x_1, \
        ") ; (x_0_g, x_1_g) = (", ws.x_0_g, ",", ws.x_1_g, ")"
    print "   (y_0, y_1) = (", ws.y_0, ",", ws.y_1, \
        ") ; (y_0_g, y_1_g) = (", ws.y_0_g, ",", ws.y_1_g, ")"
    print "   (Nx, Ny; N) = (", ws.Nx, ",", ws.Ny, ";", ws.N, ")", \
        "(Nx_g, Ny_g; N_g) = (", ws.Nx_g, ",", ws.Ny_g, ";", ws.N_g, ")"
    print "   (Nr, dr; deltat_kernel) = (", ws.Nr, ",", ws.dr, ",", \
        ws.deltat_kernel, ")"
    print "   flag_saveformat_ghostboxes =", ws.flag_saveformat_ghostboxes, ";"
    print ""
    print "   (Lx, Ly, A) = (", ws.Lx, ",", ws.Ly, ",", ws.A, ")"
    print "   (i0_g, i1_g) = (", ws.i0_g, ",", ws.i1_g, ") ; " + \
            "(j0_g, j1_g) = (", ws.j0_g, ",", ws.j1_g, ")"
    print "   (i0_f, i1_f) = (", ws.i0_f, ",", ws.i1_f, ") ; " + \
            "(j0_f, j1_f) = (", ws.j0_f, ",", ws.j1_f, ")"
    print "   (i0_i, i1_i) = (", ws.i0_i, ",", ws.i1_i, ") ; " + \
            "(j0_i, j1_i) = (", ws.j0_i, ",", ws.j1_i, ")"
    print "   flag_saveformat_text =", ws.flag_saveformat_text, \
            "; flag_saveformat_binary =", ws.flag_saveformat_binary, ";"
    print ""

###############################################################################

    
def read_array(ws, field, it=-1, cropped=1, tipo='float64', verbose=1):
    
    if it < 0:
        it = ws.Niter + it + 1
            
    if field=='ros' and ws.flag_ros_const:
        print "...WARNING (read_array): ros is constant, so loading data for it=1..."
        it = 1
    
                    
    if (ws.flag_saveformat_binary):
        
        filename = ws.path_ws + ws.filename_ws + '_' + field + '_' + \
            str(it).zfill(3) + '.' + ws.filename_extension_bin
            
        try:
            f = open(filename, 'r')
            (n, m) = np.fromfile(f, dtype=np.int32, count=2)
        except:
            print "...ERROR (read_array): '{}' not available.".format(filename)
            return None
        
        A = np.ndarray([n,m]) ;
        
        for i in range(n):
            if tipo == 'float64':
                A[i,:] = np.fromfile(f, dtype=np.float64, count=m) ;
            elif tipo == 'int32':
                A[i,:] = np.fromfile(f, dtype=np.int32, count=m) ;
                
    elif (ws.flag_saveformat_text):
        
        filename = ws.path_ws + ws.filename_ws + '_' + field + '_' + \
            str(it).zfill(3) + '.' + ws.filename_extension_txt
                
        A = []
        try:
            f = open(filename,'r')
        except:
            print "...ERROR (read_array): '{}' not available.".format(filename)
            return None
        for line in f.readlines():
            for element in line.split():
                if tipo == 'float64':
                    A.append(float(element))
                elif tipo == 'int32':
                    A.append(int(element))
        f.close()
        
        (n, m) = (int(A[0]), int(A[1]))
        A = np.asarray(A[2:]).reshape(n, m)
        
    else:
        
        print "...ERROR (read_array): '{}' not available.".format(filename)
        return None
        
    if cropped:
        A = crop_array(ws, A, verbose)
        (n, m) = np.shape(A)
        
    if verbose:
        print "...read array {}x{} ({}) from file '{}'".format(n, m, tipo, filename)
    
    return A
    

   
###############################################################################
    
def crop_array(ws, A, verbose=1):

    (n0, m0) = np.shape(A)

    (i0, i1, j0, j1) = (ws.i0_f, ws.i1_f, ws.j0_f, ws.j1_f)
    A = A[j0:j1+1, i0:i1+1]

    if verbose:
        (n, m) = np.shape(A)
        print "...array cropped: ({}x{})->({}x{}).".format(n0, m0, n, m)
            
    return A
    

        
###############################################################################
    
def read_array_gp(filename, path_ws='./', tipo='float64', verbose=1):
            
    if filename[-4:] != ".txt": filename = filename + ".txt"

    try:
        f = open(path_ws+filename,'r')
        A = [] 
        for line in f.readlines():
            for element in line.split():
                A.append(element)
        f.close()
    except:
        print "...ERROR (read_array_gp): '{}' not available.".format(path_ws+filename)
        return None

    n = (int)(np.round(np.sqrt(np.size(A))))
    
    A = np.asarray(A).reshape(n, n).astype(tipo) # nrow, ncol    
    if verbose:
        print "...read array {}x{} ({}) from file '{}{}'.\n".format(n, n, \
        path_ws, filename)
    
    return A
    

##############################################################################


def read_t(ws, it=-1, verbose=1):
    
    if it < 0: 
        it = ws.Niter + it + 1
    
    if (ws.flag_saveformat_binary):
        
        filename = ws.path_ws + ws.filename_ws + "_" + ws.filename_data + \
            "_" + str(it).zfill(3)  + '.' + ws.filename_extension_bin
            
        try:
            f = open(filename, 'r')
        except:
            print "...ERROR (read_t): '{}' not available.".format(filename)
            return None
        Nit = np.fromfile(f, dtype=np.int32, count=1)
        t = np.fromfile(f, dtype=np.float64, count=1)
        Nit, t = int(Nit), float(t)
                
    elif (ws.flag_saveformat_text):

        filename = ws.path_ws + ws.filename_ws + "_" + ws.filename_data + \
            "_" + str(it).zfill(3)  + '.' + ws.filename_extension_txt
            
        try:
            f = open(filename,'r')
        except:
            print "...ERROR (read_t): '{}' not available.".format(filename)
            return None
        data = f.read().split()
        (Nit, t) = (int(data[0]), float(data[1]))
                
    else:
        
        print "\n...Oh boy, I ain't got no file :(\n\n"
        return None

    f.close()
    
    if verbose:
        print "...Nit = {}, t = {}.".format(Nit, t)
    
    return (Nit, t)    
      

###############################################################################

def make_contour(ws, field, it=-1, show_frames=True, L=L1, cm=None, \
    cropped=1, grid=0, dim=None, tipo='float64', verbose=1):
        
    if field in [ 'tburn_phi', 'tburn_phieff', 'tburn_psi' ] and it != -1:
        it = -1
        print "...WARNING (make_contour): field '{}' only available for it=ws.Niter".format(field)
        
    if cropped:
        xx, yy = ws.xxc, ws.yyc #(xx, yy) = read_grid(ws, cropped, verbose=0) ; 
    else:
        xx, yy = ws.xx, ws.yy
        
    files = []
    os.system("rm -f frame_*.png")
        
    if isinstance(it, int):
        
        if it < 0:
            frames = [ ws.Niter + it + 1 ]
        
        elif it == 0:
            frames = ws.Niters
                
        else:
            frames = range(it, ws.Niter, it)
            frames.insert(0, 1)
            frames.append(ws.Niter)
 
    if dim is None:
        Lx, Ly = 5*xx[-1]/yy[-1], 5
    else:
        Lx, Ly = dim[0], dim[1]
        
    if not show_frames:
        plt.ioff()
    
    j = 0
    Nframes = len(frames)
    
    for iframe in frames:

        j += 1
        fig = plt.figure(figsize=(Lx, Ly))
        ax = fig.add_subplot(111)
        
        A = read_array(ws, field, iframe, cropped, tipo, verbose=0) ; 
        if A is None:
            j -= 1
            continue

        (Nit, t) = read_t(ws, iframe, verbose=0) ;

        ax.cla()
        if L==[]:
            cs = ax.contour(xx, yy, A, colors=cm)
        else:
            cs = ax.contour(xx, yy, A, levels=L, colors=cm)
        plt.title('[' + str(j) + '] it: ' + str(Nit) + ', t: ' + str(t))
        #plt.title('[{:03}] it: {:f}, t: {:f}'.format(j, Nit, t))
        plt.grid('on')
        if ws.flag_fire_obstacle_1:
            x0, x1 = max(ws.x0_fire_obstacle_1,xx[0]), min(ws.x1_fire_obstacle_1,xx[-1])
            y0, y1 = max(ws.y0_fire_obstacle_1,yy[0]), min(ws.y1_fire_obstacle_1,yy[-1])            
            #plt.add_patch(Rectangle((x0,y0),x1-x0,y1-y0, color=col_grigio30))
            ax.axhspan(y0, y1, x0/xx[-1], x1/xx[-1], facecolor=col_grigio30, edgecolor=col_grigio30 )
        if ws.flag_fire_obstacle_2:
            x0, x1 = max(ws.x0_fire_obstacle_2,xx[0]), min(ws.x1_fire_obstacle_2,xx[-1])
            y0, y1 = max(ws.y0_fire_obstacle_2,yy[0]), min(ws.y1_fire_obstacle_2,yy[-1])  
            ax.axhspan(y0, y1, x0/xx[-1], x1/xx[-1], facecolor=col_grigio30, edgecolor=col_grigio30 )
            #plt.add_patch(Rectangle((x0,y0),x1-x0,y1-y0, color=col_grigio30))
        if grid:
            # NB: disegna i bordi delle celle! (I valori sono calcolati nei centroidi)
            [ plt.plot([xn,xn],[ws.yyn[0],ws.yyn[-1]],'k') for xn in ws.xxn ]
            [ plt.plot([ws.xxn[0],ws.xxn[-1]],[yn,yn],'k') for yn in ws.yyn ]
            #[ plot([ws.yyn[0],ws.yyn[-1]],[xn,xn],'k') for yn in ws.yyn ]
        ax.set_aspect('equal')
        #axis('equal')
        ax.axis([xx[0], xx[-1], yy[0], yy[-1]])
        
        fname = 'frame_' + field + '_' + str(j).zfill(3) + '.png'
        fig.savefig(fname)
        files.append(fname)
        
        if show_frames:
            plt.draw()
            plt.pause(0.5)
        else:
            plt.close('all')
        
        if verbose:
            print "...contour '{}' created [{}/{}]".format(fname, j, Nframes)
    
    if show_frames:
        plt.ion()
    else:
        plt.close('all')
    
    return (files, fig, cs)
    
    
###############################################################################


def make_contour_movie(ws, field, frames=0, show_frames=False, L=L1, cm=None, \
    cropped=1, grid=0, dim=None, tipo='float64', verbose=1, mname='animation'):
        
    if mname == 'animation':
        os.system("rm -f animation.mpg")
        
    mname = mname +".mpg"
    
    if verbose:
        print "Making movie '"+mname+"'. It could take a while..."
        
    (files, fig, cs) = make_contour(ws, field, frames, show_frames, L, cm, \
        cropped, grid, dim, tipo, verbose) ;
        
    os.system("mencoder 'mf://frame_*.png' -quiet -mf type=png:fps=10 -ovc "+\
        "lavc -lavcopts vcodec=wmv2 -oac copy -o "+mname)

    if verbose:
        print "Movie '{}' has been created.".format(mname)
    
    return None
    
    
    
###############################################################################


def make_slice_ppp(ws, flag_xy=0, nm=-1, it=-1, flag_average=0):
    
    if it < 0: it = ws.Niter
        
    x = ws.xx
    
    phi = read_array(ws, 'phi', it) ;
    phieff = read_array(ws, 'phieff', it) ;
    psi = read_array(ws, 'psi', it) ;
    
    (ncol, nrow) = phi.shape 
    even_col = 1 - np.mod(ncol,2)
    even_row = 1 - np.mod(nrow,2)
    
    if flag_xy==0:
        
        m = 0
        
        if nm <= 0: 
            n = int(ncol/2) + np.mod(ncol,2) - 1
        else:
            n = nm
            flag_average = 0
            
        if flag_average and even_row:
            strn, strm = str(n+0.5), ':'
            Phi = (2.0-phi[n,:]-phi[n+1,:])/2.0
            Phieff = (phieff[n,:]+phieff[n+1,:])/2.0
            Psi = (psi[n,:]+psi[n+1,:])/2.0
        else:
            strn, strm = str(n), ':'
            Phi = 1.0-phi[n,:]
            Phieff = phieff[n,:]
            Psi = psi[n,:] ;
            
    else:
        
        n = 0
        
        if nm <= 0: 
            m = int(nrow/2) + np.mod(nrow,2) - 1
        else:
            m = nm
            flag_average = 0
            
        if flag_average and even_col:
            strn, strm = ':', str(m+0.5)
            Phi = (2.0-phi[:,m]-phi[:,m+1])/2.0
            Phieff = (phieff[:,m]+phieff[:,m+1])/2.0
            Psi = (psi[:,m]+psi[:,m+1])/2.0
        else:
            strn, strm = ':', str(m)
            Phi = 1.0-phi[:,m]
            Phieff = phieff[:,m] 
            Psi = psi[:,m] 
    
    f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
    ax1.plot(x, Phi, '.-k')
    ax1.set_title('(it:'+str(it)+'): 1-phi / phieff / psi ['+strn+','+strm+']')
    ax1.axis([0, 10000, 0, 1.1])
    ax2.plot(x, Phieff,'.-b')
    ax3.plot(x, Psi, '.-r')
    # Fine-tune figure; make subplots close to each other and hide x ticks for
    # all but bottom plot.
    f.subplots_adjust(hspace=.1)
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    
    return None
    
####################################################################
    
def readfromfile(fname, flag_array=False):
    
    f = open(fname,'r')
    (beta, N) = f.readline().split()
    beta, N = float(beta), int(N)
    rr, yy = [], []
    for i in range(N):
        (r, y) = f.readline().split()
        rr.append(float(r))
        yy.append(float(y))
        
    f.close()
    if flag_array:
        rr, yy = np.asarray(rr), np.asarray(yy)
    
    return (rr, yy)
    
def get_M_data_radial(beta, flag_array=True):
    
    if beta == 0.25:
        filename = './adlib/M_data_radial_beta_1_4.txt'
    elif beta == 0.5:
        filename = './adlib/M_data_radial_beta_1_2.txt'
    elif beta == 0.75:
        filename = './adlib/M_data_radial_beta_3_4.txt'
        
    (rr, yy) = readfromfile(filename, flag_array)
    print "read data from file '" + filename + "'."
    
    return (rr, yy)
    
def makefig_M_data(ifig=None):
    (rr1, yy1) = get_M_data_radial(0.25)
    (rr2, yy2) = get_M_data_radial(0.5)
    (rr3, yy3) = get_M_data_radial(0.75)
    
    plt.figure(num=ifig)
    plt.plot(rr1, yy1, '.-b', rr2, yy2, '.-k', rr3, yy3, '.-r')
    plt.axis([0, 5, 0, 3.5])
    plt.grid(True)
    
    return None