import numpy as np
import matplotlib.pyplot as plt
import scipy
import sys
import field_calc

# Reads a binary file and plots a cut in the x-y plane.
# Execute in ipython with '%run plot_bindata.py'

def readslice(inputfilename,nx,ny,nfile):
   f = open(inputfilename,'rb')         #Opens file to read (and also accept binary output I think)
   f.seek(4)                            #Fortran files begin with integer telling you how many bytes are recorded in the write. We skip this integer and go straight to the data.
   field = np.fromfile(f,dtype = 'd',count=nx*ny/nfile)         #gets the data in a 1x(nx*ny/nfile) vector, later to be transformed into 2D array, once stitched together
   f.close()
   return field


# Path to the binary data
runname = raw_input("Folder name: ")
path = '../'+runname+'/outs/'

if ('kf24' in runname) or ('1024' in runname):
	num_files = 32 #input('Enter number of output files (abc.xyz.out, give highest abc+1): ')
	reso = 1024#input("Resolution? :")
elif ('kf6' in runname) or ('256' in runname):
	num_files = 16
	reso = 256
else:
	num_files=16
	reso = 512

# Spatial resolution
NX = reso
NY = reso
shape = (NX,NY)

tf = np.loadtxt('../'+runname+'/run/time_field.txt')
print("Last output: %s" % int(tf[-1][0]))

outnum = raw_input("out num? ") #sys.argv[1]
outnum ="{:0>3s}".format(outnum)

#otypes = ['ps','ww','vz','vx','vy']
otypes = ['ww','vy','vz']
legends={'ps':r'$\psi$','ww':r'$\omega_z$','vz':r'$v_z$','vy':r'$v_y$','vx':r'$v_x$'}

# Reads binary files
#psi = np.fromfile(path+'ps.'+outnum+'.out',dtype=np.float32).reshape(shape,order='F')
outs  = dict([])
datbars = dict([])
for otype in otypes:
	out = field_calc.field_calc(runname,otype,outnum,reso=reso,num_files=num_files)
	outs[otype] = out

	omax = np.max(out)
	omin = np.min(out)
	datbars[otype] = np.max([abs(omax),abs(omin)])
	#datbars[otype] = 250

# Show a horizontal cut of the field in the middle of the box
for otype in otypes:
	datbar = datbars[otype]
	plt.figure()
	#plt.imshow(outs[otype],vmin = -datbar,vmax = datbar,cmap='bwr')
	plt.imshow(outs[otype],vmin = -datbar,vmax = datbar,cmap='bwr')
	plt.title(runname+' '+legends[otype]+' out = '+str(outnum))
	#plt.title(runname+' '+otype+' out = '+str(outnum))
	plt.colorbar()
	plt.xlabel('x')
	plt.ylabel('y')
	plt.tight_layout()
	plt.savefig('./figs/'+runname+'_'+otype+'_'+str(outnum)+'.png',bbox_inches='tight',dpi=300)
	
plt.show()
