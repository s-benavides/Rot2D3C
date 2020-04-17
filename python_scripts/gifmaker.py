##############  REMEMBER TO RUN THIS IN THE D6S1Q10/vis_(outtype)/ folder!  ###############

### Then to make gif type: " convert   -delay 10   -loop 0   FlowD6S1Q100_ww_*.png   D6S1Q100.gif" into command line

import scipy
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os, sys
from pylab import *

def readslice(inputfilename,nx,ny,nfile):
   f = open(inputfilename,'rb')         #Opens file to read (and also accept binary output I think)
   f.seek(4)                            #Fortran files begin with integer telling you how many bytes are recorded in the write. We skip this integer and go straight to the data.
   field = np.fromfile(f,dtype = 'd',count=nx*ny/nfile)         #gets the data in a 1x(nx*ny/nfile) vector, later to be transformed into 2D array, once stitched together
   f.close()
   return field

num_files = input('Enter number of slices (abc.xyz.out, give highest abc+1): ')
outnum_strt = input('Enter first output number (abc.xyz.out, give xyz_start): ')
outnum_end = input('Enter last output number (abc.xyz.out, give xyz_final): ')
reso = input('Resolution of simulation? (512, 1024, etc.): ')
Ntemp = raw_input('Enter N: ')
Ftemp = raw_input('Enter F: ')
Ktemp = raw_input('Enter K (dont include _vis part): ')
outtype = raw_input('Enter flow file (ps, ww, E2D): ')

outnum_strt = np.int(outnum_strt)
outnum_end = np.int(outnum_end)

tot_num = outnum_end - outnum_strt

if outtype in ['E2D']:
        for j in range(tot_num):

                data1 = []
                
                for i in range(num_files):
                        data1.append(readslice('../cond_ps.'+str("%03d" % i)+'.'+str("%03d" % (outnum_strt + j))+'.out',reso,reso,num_files)) #adds all the 1x(nx*ny/nfile) to the total list

                final1 = np.reshape(data1,(reso,reso)) #Finally reshapes into a reso x reso array.

		########### Below is for plotting energy of 2D velocities:

		dx,dy = gradient(final1)
		EnP = [[dx[k][l]**2 + dy[k][l]**2 for k in range(shape(dx)[0])] for l in range(shape(dx)[0])]
		
                if j==0:
                        datmin = min(min(EnP))
                        datmax = max(max(EnP))

		figure(1)
		im1 = plt.imshow(EnP, cmap=cm.hot)
		if j==0:
			cbar1 = plt.colorbar(im1)
		plt.title('N'+Ntemp+'F'+Ftemp+'K'+Ktemp+' 2D Energy')
		plt.savefig('FlowN'+Ntemp+'F'+Ftemp+'K'+Ktemp+'_E2D_'+str("%03d" % (j+1))+'.png')

else:
	for j in range(tot_num):

        	data1 = []

       		for i in range(num_files):
               		 data1.append(readslice('../cond_'+outtype+'.'+str("%03d" % i)+'.'+str("%03d" % (outnum_strt + j))+'.out',reso,reso,num_files)) #adds all the 1x(nx*ny/nfile) to the total list

		final1 = np.reshape(data1,(reso,reso)) #Finally reshapes into a reso x reso array.
		
		if j==0:
			datmin = min(min(final1))
		#	print(datmin)
			datmax = max(max(final1))
		#	print(datmax)	
		#datmin = -50
		#datmax = 50
              
        	figure(1)
     		im1 = plt.imshow(final1, cmap=cm.hot,vmin=datmin,vmax=datmax)
        	if j==0:
                	cbar1 = plt.colorbar(im1)
        	plt.title('N'+Ntemp+'F'+Ftemp+'K'+Ktemp+'_'+outtype+'_'+str("%03d" % (j+1)))
      		plt.savefig('FlowN'+Ntemp+'F'+Ftemp+'K'+Ktemp+'_'+outtype+'_'+str("%03d" % (j+1))+'.png')
