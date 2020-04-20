import numpy as np
import glob as glob
import string as string
import matplotlib.pyplot as plt

# Reads binary files in a directory and does
# postprocessing (e.g., to visualize later with
# VAPOR). Note that for big runs, GHOST can do some
# automatic postprocessing (e.g., compute vorticity)
# at run time.
# Execute in ipython with '%run postprocess.py'

def readslice(inputfilename,nx,ny,nfile):
   f = open(inputfilename,'rb')         #Opens file to read (and also accept binary output I think)
   f.seek(4)                            #Fortran files begin with integer telling you how many bytes are recorded in the write. We skip this integer and go straight to the data.
   field = np.fromfile(f,dtype = 'd',count=nx*ny/nfile)         #gets the data in a 1x(nx*ny/nfile) vector, later to be transformed into 2D array, once stitched together
   f.close()
   return field


def field_calc(run,otype,outnum,reso=512,num_files=16):
	# Path to the binary data
	path = '../'+run+'/outs/'
	
	# Box size
	Lx = 2*np.pi
	Ly = 2*np.pi
	
	# Spatial resolution
	NX = reso
	NY = reso
	dx = Lx/NX
	dy = Ly/NY
	shape = (NX,NY)

	data = []
        for i in range(num_files):
                data.append(readslice(path+'ps.'+str("%03d" % i)+'.'+outnum+'.out',reso,reso,num_files)) #adds all the 1x(nx*ny/nfile) to the total list

        ps = np.reshape(data,(reso,reso)) #Finally reshapes into a reso x reso array.

	# Reads binary files, computes vertical vorticity
	# using centered finite differences, and saves in 
	# a new binary file named 'wz.NNNN.out'
	if 'vy' in otype:
		adv = np.roll(ps,-1,axis=1)
		ret = np.roll(ps,1,axis=1)
		output = -(adv-ret)/(2*dx) # -dpsi/dx
	elif 'vx' in otype:
                adv = np.roll(ps,-1,axis=0)
                ret = np.roll(ps,1,axis=0)
                output = (adv-ret)/(2*dy) # dpsi/dy
	elif 'vz' in otype:
	        data = []
	        for i in range(num_files):
        	        data.append(readslice(path+'vz.'+str("%03d" % i)+'.'+outnum+'.out',reso,reso,num_files)) #adds all the 1x(nx*ny/nfile) to the total list

        	output = np.reshape(data,(reso,reso)) #Finally reshapes into a reso x reso array.
        elif 'ww' in otype:
                data = []
                for i in range(num_files):
                        data.append(readslice(path+'ww.'+str("%03d" % i)+'.'+outnum+'.out',reso,reso,num_files)) #adds all the 1x(nx*ny/nfile) to the total list

                output = np.reshape(data,(reso,reso)) #Finally reshapes into a reso x reso array.
	else:
		output = ps
	
	return output
