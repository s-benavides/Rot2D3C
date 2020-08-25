import matplotlib
matplotlib.use('Agg')
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import glob as glob
import sys
import imageio
import field_calc
#imageio.plugins.ffmpeg.download()

plt.style.use('default')
plt.rcParams.update({'font.size': 20})
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams["font.serif"] = "Times New Roman"

# Reads a binary file and plots a cut in the x-y plane.
# Execute in ipython with '%run plot_bindata.py'

# Path to the binary data
runname = raw_input("Folder name: ")
path = '../'+runname+'/outs/'

# Out type (e.g. 'wz','ax','spec2D_yavg',...)
#otype = 'kspec2D_yavg'
otype = 'ww'
#otype = 'vy'
#otype = 'vz'
print("Making movie of %s for run %s" % (otype,runname))

# Spatial resolution
if ('kf24' in runname) or ('1024' in runname):
        num_files = 32 #input('Enter number of output files (abc.xyz.out, give highest abc+1): ')
        N = 1024#input("Resolution? :")
elif ('kf6' in runname) or ('256' in runname):
        num_files = 16
        N = 256
else:
        num_files=16
        N = 512

shape = (int(N),int(N))

# File path
filelist = sorted(glob.glob(path+'ps.001.*.out'))
nfiles = np.size(filelist)
print("nfiles = %s" % nfiles)
 
frames_per_second = 5
writer = imageio.get_writer('./movies/'+runname+'_'+otype+'.mp4', mode="I", codec='h264',bitrate=1800000,fps = frames_per_second,quality=10.0)#output_params=['-s','1500x1500'])
#writer = imageio.get_writer('./movies/'+runname+'_'+otype+'.wmv', codec='msmpeg4',quality=10.0, fps = frames_per_second)

legend_otype = {'ww':r'$\omega$','vy':r'$v_y$','vz':r'$v_z$'}

data = pickle.load(open(str('./rundat/Rot2D3C_DataAvg_2020-06-29.p'),'rb'))
omegaz = data[runname]['omega']
kf = data[runname]['kf']
inj, dinj,_ = data[runname]['inj']  # we don't force the magnetic field
uf = (5./4.)*(inj/kf)**(1./3.)
OoRo= (2*omegaz)/(kf*uf)

# Reads binary files
for ii,ofile in enumerate(filelist):
	# Reads binary files
	ind = ofile.split(path+'ps.001.')[1]
	ind = ind.split('.out')[0]
	print("Working on snapshot %s, out num %s" % (ii+1,ind))
        out = field_calc.field_calc(runname,otype,ind,reso=N,num_files=num_files)

	datmin = np.amin(out)
	datmax = np.amax(out)
	datbar = max(abs(datmin),abs(datmax))

	fig = plt.figure(1,figsize=(12,10))
	im = plt.imshow(out,cmap=cm.bwr,vmin=-datbar,vmax=datbar)
	cbar = plt.colorbar(im)
	plt.xlabel(r'$x$')
	plt.ylabel(r'$y$')
	plt.title("Vorticity "+legend_otype[otype]+r", $Ro^{-1} =$ %.2f, $k_f = $ %s, out # %s" % (OoRo,kf,ind))
	plt.tight_layout()
	fig.canvas.draw()
	img = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
	img = img.reshape(fig.canvas.get_width_height()[::-1] + (3,))
	writer.append_data(img)
	plt.close()

print("done making %s for run %s" % (otype,runname))
