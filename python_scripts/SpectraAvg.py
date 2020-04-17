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
avglist = sorted(glob.glob('rundat/*T90*.txt'))

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
        [start,start_fl] = np.loadtxt('rundat/AvgTime'+run+'.txt')

        start = int(start)
        start_fl = int(start_fl)

	rinfo  = np.genfromtxt('../'+run+'/run/parameter.inp',comments='!',skip_footer=136,skip_header=15,converters={2:  lambda val: float(val.translate(rule))},usecols=2)

        rand = rinfo[5]	

	inds = []
	files = sorted(glob.glob(path+'kspectrum.*.txt'))
        for file in files:
                num = file.split(path+'kspectrum.')[1]
                num = num.split('.txt')[0]
                if (int(num)>=start_fl):
                        inds.append(num)


        numfiles = int(len(inds))
        print('%s files to average' % numfiles)


        # Reads flux files
        if rand!=0:
        	injtot = 0.5*injtot

	fields = ['kspectrum','mspectrum','kspecperp','mspecperp','kspecpara','mspecpara']

	# Averaging
	for field in fields:
        	for ii,ind in enumerate(inds):
	        	# Load file names
                	flux = sorted(glob.glob(path+field+'.'+str(ind)+'.txt'))[0]
			# Load and Average File Names:
                	if ii==0:
                		flux_avg = np.loadtxt(flux)[:,1]/float(numfiles)
			else:
                		flux_avg += np.loadtxt(flux)[:,1]/float(numfiles)
		
		Data_E[field] = flux_avg
	
	Data_E['ks'] = np.loadtxt(flux)[:,0]
        Data[run] = Data_E

# Saving data:
print "Saving Data"
name = 'rundat/SpectraAvg_'+str(date.today())+'.p'
pickle.dump(Data, open(str(name), 'wb'))
