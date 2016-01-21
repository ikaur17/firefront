# -*- coding: utf-8 -*-

import forefire, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.cm as cm
import subprocess
import pylab
import matplotlib.mlab as ml
import scipy.interpolate

def getLocationFromLine(line):

    llv = line.split("loc=(")
    if len(llv) < 2: 
        return None
    llr = llv[1].split(",");
    if len(llr) < 3: 
        return None
    return (float(llr[0]),float(llr[1]))


def dist(a,b):
    return np.sqrt(np.power(a[0]-b[0],2)+np.power(a[1]-b[1],2))


def printToPathe(linePrinted):
    fronts = linePrinted.split("FireFront")
    pathes = []
    for front in fronts[1:]:
        nodes =front.split("FireNode")[1:]
        if len(nodes) > 0: 
            Path = mpath.Path
           
            codes = []
            verts = []
            lastNode = getLocationFromLine(nodes[0])
            firstNode = getLocationFromLine(nodes[0])
            
 
            codes.append(Path.MOVETO)
            verts.append(firstNode)      

            for node in nodes[:]:
                newNode = getLocationFromLine(node)
                codes.append(Path.LINETO)
                verts.append(newNode)         
                lastNode = newNode
                
    
            codes.append(Path.LINETO)
            verts.append(firstNode)          
           
            pathes.append(mpath.Path(verts, codes))

    return pathes;

def addToPathe(xy):
	pathes = []
	Path = mpath.Path
	codes = []
	verts = []
	
	lastNode = xy[0,:]
	firstNode = xy[0,:]

	print lastNode, firstNode

	codes.append(Path.MOVETO)
        verts.append(firstNode)
#	print 'length of xy',len(xy)
	for  i in range(0,len(xy)):
			newNode = xy[i,:]
 			codes.append(Path.LINETO)
                	verts.append(newNode)
                	lastNode = newNode

	codes.append(Path.LINETO)
        verts.append(firstNode)

        pathes.append(mpath.Path(verts, codes))

	return pathes

ff = forefire.PLibForeFire()
    
###############################################################################        

sizeX = 250
sizeY = 250
scale = 20
ff.setString("fuelsTableFile","fuels.ff")
ff.setString("ForeFireDataDirectory","test")

ff.setDouble("spatialIncrement",2)
#ff.setDouble("spatialIncrement",4)
ff.setDouble("minimalPropagativeFrontDepth",1)
ff.setInt("atmoNX",sizeX)
ff.setInt("atmoNY",sizeY)

ff.setDouble("initialFrontDepth",1)
ff.setDouble("perimeterResolution",8)
#ff.setDouble("perimeterResolution",18)
ff.setDouble("relax",.2)
ff.setDouble("smoothing",5)
ff.setDouble("z0",0.)
ff.setDouble("windU",0.)
ff.setDouble("windV",3)
ff.setInt("defaultFuelType",9)
ff.setInt("defaultHeatType",0)
ff.setDouble("nominalHeatFlux",40000)
ff.setDouble("burningDuration",0)
ff.setDouble("maxFrontDepth",12)
#ff.setDouble("maxFrontDepth",20)
#ff.setDouble("minSpeed",0.000000001)
ff.setDouble("minSpeed",0.01)


ff.execute("FireDomain[sw=(0.,0.,0.);ne=(5000,5000,0.);t=0.]")

print "resolution of bmap is ", ff.getString("bmapResolution")

ff.addLayer("data","altitude","z0")
ff.addLayer("data","windU","windU")
ff.addLayer("data","windV","windV")
ff.addLayer("BRatio","BRatio","BRatio")
ff.addLayer("flux","heatFluxBasic","defaultHeatType")
ff.addLayer("propagation","BalbiUnsteady","propagationModel")
fuelmap = np.zeros((sizeX,sizeY,1), dtype=np.int32)
#fuelmap = np.zeros((8000,8000,1), dtype=np.int32)
fuelmap[:,:,:] = 1
fuelmap[:,155:158,:] = 8
#fuelmap[:,180:184,:] = 8
#fuelmap[:,155:157,:] = 8
#fuelmap[70:80,:,:] = 8
ff.addIndexLayer("table","fuel",0 , 0, 0, sizeX*scale, sizeY*scale, 0, fuelmap)

print(np.shape(ff.getDoubleArray("fuel")))

path_ff = os.getcwd()
path_ls = os.path.join(os.getcwd(), "LSfire")
path_data = os.path.join(os.getcwd(), "data/FF_LS/")
os.chdir(path_data)
subprocess.call('rm *.txt' , shell = True)
os.chdir(path_ff)

###############################################################################

def run(ff, step, Niter, t=0):
    
    print "*** running simulation: step={}, Niter={} ***".format(step, Niter)
    ff.execute("\tFireFront[t="+str(np.float(t))+"]")
    
    x, y, pathes, path1, fline, theta = [], [], [], [], [],[]
    path_old =[]
    xfinal, yfinal, tfinal = [], [], []
    icontour = 0
    for i in range(0, Niter+1):    
        print "goTo[t=%f]"%(i*step)
        ff.execute("goTo[t=%f]"%(i*step))
        pathes += printToPathe( ff.execute("print[]"))
        path1 = printToPathe( ff.execute("print[]"))
	index, x, y, fline = [], [], [], [] ; ii = 0
	x1, y1 = [], []
	xi, yi, zi = [], [], []
	xlast , ylast = [], []
    	#for path in pathes:
	#print 'len of path1',len(path1)
    	for path in path1:
	    ii = ii+1
            v=path.vertices
            x.extend(v[:,0])
            y.extend(v[:,1])
	    index.extend(v[:,0]/v[:,0]*ii)

#	    if i*step % 1200.0 == 0   :
	    if i*step % 1221.0 == 0   :
		if len(path1) == 1:
	    		xfinal.extend(v[:,0])
	    		yfinal.extend(v[:,1])
            		tfinal.extend(v[:,0]/v[:,0]*(i*step))
		elif ii == len(path1):
				xfinal.extend(v[:,0])
                        	yfinal.extend(v[:,1])
                        	tfinal.extend(v[:,0]/v[:,0]*(i*step))
		elif ii == 2:
			if max(v[:,1]) >= 3160.0:
				xfinal.extend(v[:,0])
                       		yfinal.extend(v[:,1])
                       		tfinal.extend(v[:,0]/v[:,0]*(i*step))

	np.savetxt("%s/xarray.txt" % path_data,x,fmt='%10.5f')
        np.savetxt("%s/yarray.txt" % path_data,y,fmt='%10.5f')
        np.savetxt("%s/xyarray.txt" % path_data,np.c_[x,y,index],fmt='%f')
        
	print len(x),len(y),len(fline)
        if i > -1 :
		print "Including turbulence from LS at time =", i*step
		time = str(i*step)
		os.chdir(os.path.abspath(path_ls))
        	subprocess.call('./FF_LS.sh  "%s"'%time , shell = True)
		os.chdir(path_ff)
		num_lines = sum(1 for line in open(path_data + 'new_ignition.txt'))
		if num_lines > 10:
			x1 , y1 = [], []
			xx , yy, xnew, ynew = [], [], [], []
			x1, y1 = np.loadtxt(path_data + "new_ignition.txt", usecols = (0,1), unpack = True)


			theta = np.arctan2(y1-2500.0,x1-2500.0)

			xx = zip(*sorted(zip(theta,x1), reverse = True))[1]
			#xx = zip(*sorted(zip(theta,x1)))[1]
			yy = zip(*sorted(zip(theta,y1), reverse = True))[1]
			#yy = zip(*sorted(zip(theta,y1)))[1]
			
			
			ff.execute("\tFireFront[t="+str((i+1)*step)+"]")
                        for j in range (0,len(xx)):
                                     	theta1 = np.arctan2(yy[j]-2500.0,xx[j]-2500.0)
                                	vx = np.cos(theta1)
                		        vy = np.sin(theta1)
		                        
					call_FireNode(loc=(xx[j],yy[j],0), vel=(vx,vy,0), t=(i+1)*step)

			


			

        simdata = ff, x, y, pathes, path1, fline
    xp = np.asarray(xfinal)
    yp = np.asarray(yfinal)
    zp = np.asarray(tfinal)
    np.savetxt(path_data + 'final_points.txt',np.c_[xfinal,yfinal,tfinal],fmt='%f')
    fig, ax = plt.subplots()
#    plt.colorbar()
    tab = np.transpose(ff.getDoubleArray("fuel"))[0]
    CS = ax.imshow(tab, extent=[0,5000,0,5000], aspect='auto', origin = 'lower', interpolation='nearest', cmap='gray_r')
    ax.scatter(xfinal, yfinal, c=tfinal, marker = 'o', s=3,edgecolors='none', cmap='jet')
    plt.show()


    return simdata


def call_FireNode(loc, vel=(0,0,0), t=0):

    x, y, z = loc
    u, v, w = vel
    def S(a): return str(np.float(a)) 
    s = "\t\tFireNode[loc=("+S(x)+","+S(y)+","+S(z)+");"+\
        "vel=("+S(u)+","+S(v)+","+S(w)+");t="+S(t)+"]"
    ff.execute(s)
    
    return None


    
def init_circle(ff, center=(0.0,0.0), r=1.0, n=100):

    xc, yc = center
    dtheta = 2 * np.pi / (n-1)
    for i in range(n-1):
        x, y = xc+np.cos(-dtheta*i)*r, yc+np.sin(-dtheta*i)*r
	theta1 = np.arctan2(x-xc,y-yc)
	vx = np.cos(theta1)
	vy = np.sin(theta1)
        call_FireNode(loc=(x,y,0), vel=(0,0,0), t=0)
    print "*** initialized circle: C=({},{}), r={} ***".format(xc,yc,r)

    return None
    
def init_square(ff, center=(3000,3000), l=200):

    x0, y0 = center
    call_FireNode(loc=(x0-l/2.,y0-l/2.,0))
    call_FireNode(loc=(x0-l/2.,y0+l/2.,0))
    call_FireNode(loc=(x0+l/2.,y0+l/2.,0))
    call_FireNode(loc=(x0+l/2.,y0-l/2.,0))
    print "*** initialized square: C=({},{}), l={} ***".format(x0,y0,l)
    
    return None
    


def make_contour(simdata, extent=[0, 5000, 0, 5000], fig_name=None):
    
    ff, x, y, pathes, path1, fline = simdata
    #shape = (200, 200)
    #print pathes
    fig, ax = plt.subplots()
    tab = np.transpose(ff.getDoubleArray("fuel"))[0]
    #CS = ax.imshow(tab, cmap=plt.cm.gray,extent=[0,6000,0,6000], aspect='auto', interpolation='nearest')
    cmap = plt.cm.gray
    cmap.set_bad('white')
    CS = ax.imshow(tab, extent=[0,5000,0,5000], aspect='auto', origin = 'lower', interpolation='nearest')
    #CS = ax.imshow(tab, cmap=plt.cm.gray,interpolation='nearest')
    cbar = plt.colorbar(CS)
    #cbar.ax.set_ylabel('v')
    
    for path in pathes:
        patch = mpatches.PathPatch(path,edgecolor='red',facecolor='none',  alpha=1)
        ax.add_patch(patch)
        
    ax.grid()
    ax.axis('equal')
    plt.grid(b=True, which='major', color='w', linestyle='-')
    plt.show()

    
    #fig_name = 'fuel3.png'
    if fig_name is not None:
        fig.savefig(figname ,bbox_inches='tight')
        
    # extract the contour values from the vertices and codes matrix


def save_data(simdata):
   """save data to a text file"""
   ff, x , y, pathes, path1, fline = simdata
   np.savetxt('xarray_18.txt',x,fmt='%10.5f')
   np.savetxt('yarray_18.txt',y,fmt='%10.5f')
   np.savetxt('fhistory_18.txt',fline,fmt='%10.5f')
    
   return None
    
###############################################################################

if __name__ == '__main__':
    
    init_circle(ff, center=(2500,2500), r=300, n=100)
#    init_square(ff, center=(2500,2500), l=600)

#    S = run(ff, 200, 25 )
    S = run(ff, 111, 60 ) ### the time step should be equal to the CFL criteria 

#    make_contour(S, extent=[0, sizeX*scale, 0, sizeY*scale])

    save_data(S)
