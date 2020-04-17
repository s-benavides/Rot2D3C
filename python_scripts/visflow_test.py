import scipy
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import os, sys
from pylab import *

def readslice(inputfilename,nx,ny,nfile):
   f = open(inputfilename,'rb')		#Opens file to read (and also accept binary output I think)
   f.seek(4)				#Fortran files begin with integer telling you how many bytes are recorded in the write. We skip this integer and go straight to the data.
   field = np.fromfile(f,dtype = 'd',count=nx*ny/nfile)  	#gets the data in a 1x(nx*ny/nfile) vector, later to be transformed into 2D array, once stitched together
   f.close()
   return field

num_files = 16 #input('Enter number of output files (abc.xyz.out, give highest abc+1): ')
outnum = input('Enter output number (abc.xyz.out, give xyz): ')
reso = 512 #input('Resolution of simulation? (512, 1024, etc.): ')
#Ntemp = raw_input('Enter N: ')
#Ftemp = raw_input('Enter F: ')
#Ktemp = raw_input('Enter K: ')
#qname = raw_input('Enter Q: ')
#outtype = raw_input('Enter flow file (vx, vy, vz, ps, ww): ')
otypes = ['ww','phi']
datbars = dict([])
for otype in otypes:
	data =[]
	for i in range(num_files):
		data.append(readslice('../outs/'+otype+'.'+str("%03d" % i)+'.'+str("%03d" % outnum)+'.out',reso,reso,num_files)) #adds all the 1x(nx*ny/nfile) to the total list
	
	#print(shape(data1))
	out = np.reshape(data,(reso,reso)) #Finally reshapes into a reso x reso array.
	
	out = abs(out)

	omax = np.max(out)
	omin = np.min(out)
	datbars[otype] = np.max([abs(omax),abs(omin)])
	
	
	figure()
	im = plt.imshow(out,vmin = -datbars[otype],vmax=datbars[otype], cmap=cm.bwr)
	cbar1 = plt.colorbar(im)
	plt.title(otype)
	plt.xticks(visible=False)
	plt.yticks(visible=False)
	#plt.savefig('FlowN%sF%sK%s_ps_out%s.png' % (Ntemp,Ftemp,Ktemp,outnum))
plt.show()
