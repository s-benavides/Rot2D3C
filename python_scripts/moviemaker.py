import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import glob as glob
import sys
import imageio
import field_calc
#imageio.plugins.ffmpeg.download()

# Reads a binary file and plots a cut in the x-y plane.
# Execute in ipython with '%run plot_bindata.py'

# Path to the binary data
runname = raw_input("Folder name: ")
path = '../'+runname+'/outs/'

# Out type (e.g. 'wz','ax','spec2D_yavg',...)
#otype = 'kspec2D_yavg'
#otype = 'ww'
otype = 'vy'
print("Making movie of %s for run %s" % (otype,runname))

# Spatial resolution
N = 512
shape = (int(N),int(N))
num_files = 16 # Number of cores used

# File path
filelist = sorted(glob.glob(path+'ps.001.*.out'))
nfiles = np.size(filelist)
print("nfiles = %s" % nfiles)
 
frames_per_second = 5
writer = imageio.get_writer('./movies/'+runname+'_'+otype+'.wmv', codec='msmpeg4',quality=10.0, fps = frames_per_second)


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

	fig = plt.figure(1)
	im = plt.imshow(out,cmap=cm.bwr,vmin=-datbar,vmax=datbar)
	cbar = plt.colorbar(im)
	plt.xlabel('x')
	plt.ylabel('y')
	plt.title("%s, run %s, out # %s" % (otype,runname,ind))
	plt.tight_layout()
	fig.canvas.draw()
	img = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
	img = img.reshape(fig.canvas.get_width_height()[::-1] + (3,))
	writer.append_data(img)
	plt.close()

print("done making %s for run %s" % (otype,runname))
