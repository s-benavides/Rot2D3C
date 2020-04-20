import numpy as np
import matplotlib.pyplot as plt
#from matplotlib import rc
import pickle
import glob as glob
import string
from datetime import date

# Making rule for reading nui,nun,mu
rule = string.maketrans('d', '0')

# Get all the AvgTimeO*B*.txt files
avglist = sorted(glob.glob('rundat/AvgTimeO*[!_c].txt'))

runnames = []
for file in avglist:
  run = file.split('rundat/AvgTime')[1]
  run = run.split('.txt')[0]
  runnames.append(run)

Data = dict([])
for i,run in enumerate(runnames):
        Data_E = dict([])
        print("Working on run %s " % run)
        path = '../'+run+'/outs/'
	
	 # Average start indices
        [start,start_fl,err_ind] = np.loadtxt('rundat/AvgTime'+run+'.txt')

        start = int(start)
        start_fl = int(start_fl)

	inds = []
	files = sorted(glob.glob(path+'spectrum.*.txt'))
        for file in files:
                num = file.split(path+'spectrum.')[1]
                num = num.split('.txt')[0]
                if (int(num)>=start_fl):
                        inds.append(num)


        numfiles = int(len(inds))
        print('%s files to average' % numfiles)


	# Averaging
        for ii,ind in enumerate(inds):
	       	# Load file names
               	flux = sorted(glob.glob(path+'spectrum.'+str(ind)+'.txt'))[0]
		# Load and Average File Names:
               	if ii==0:
               		flux_2D_avg = np.loadtxt(flux)[:,0]/float(numfiles)
               		flux_z_avg = np.loadtxt(flux)[:,1]/float(numfiles)
		else:
               		flux_2D_avg += np.loadtxt(flux)[:,0]/float(numfiles)
               		flux_z_avg += np.loadtxt(flux)[:,1]/float(numfiles)
	
	Data_E['E_2D'] = flux_2D_avg
	Data_E['E_z'] = flux_z_avg
        Data[run] = Data_E

# Saving data:
print "Saving Data"
name = 'rundat/Rot2D3C_SpectraAvg_'+str(date.today())+'.p'
pickle.dump(Data, open(str(name), 'wb'))
