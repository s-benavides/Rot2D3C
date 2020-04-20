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

	params = np.genfromtxt('../'+run+'/run/parameter.inp',comments='%',converters={0:  lambda val: float(val.translate(rule))},usecols=0)
	
	rand = params[21]
	
	# Reads balance.txt
	t,en,enz,inj,injz,diss,dissz,hdiss,hdissz,coup,uf,ufz = np.loadtxt('../'+run+'/run/energy_bal.txt',unpack=True)
	t,enst,injenst,dissenst,hdissenst = np.loadtxt('../'+run+'/run/enstrophy_bal.txt',unpack=True)
	
	if rand==3:
		injenst=injenst/2.
		injz = injz/2.
	        inj = inj/2. # RANDOM FORCING
	
	# Averaging injection
	injtot = np.mean(inj[start:]+injz[start:])
	injenst=np.mean(injenst[start:])
	
	# nu,kf
	nu = float(params[11])
	hnu = float(params[12])
	nn = float(params[13])
	mm = float(params[14])
	kdn = float(params[8])
	kup = float(params[10])
	kf = (kdn+kup)/2.
	
	# nuv, etc
	nuv = float(params[15])
	hnuv = float(params[16])
	nnv = float(params[17])
	mmv = float(params[18])


	inds = []
        files = sorted(glob.glob(path+'fluxes.*.txt'))
        for file in files:
                num = file.split(path+'fluxes.')[1]
                num = num.split('.txt')[0]
                if (int(num)>=start_fl):
                        inds.append(num)

        numfiles = int(len(inds))
        print('%s files to average' % numfiles)


	# Averaging
        for ii,ind in enumerate(inds):
	       	# Load file names
               	flux = sorted(glob.glob(path+'fluxes.'+str(ind)+'.txt'))[0]
               
		# Load and Average File Names:
		if ii==0:
               		flux_enst_avg = np.loadtxt(flux)[:,0]/float(numfiles)/injenst
               		flux_2D_avg = np.loadtxt(flux)[:,1]/float(numfiles)/injtot
               		flux_z_avg = np.loadtxt(flux)[:,2]/float(numfiles)/injtot
		else:
               		flux_enst_avg += np.loadtxt(flux)[:,0]/float(numfiles)/injenst
               		flux_2D_avg += np.loadtxt(flux)[:,1]/float(numfiles)/injtot
               		flux_z_avg += np.loadtxt(flux)[:,2]/float(numfiles)/injtot
	
	Data_E['Z'] = flux_enst_avg
	Data_E['E_2D'] = flux_2D_avg
	Data_E['E_z'] = flux_z_avg
        Data[run] = Data_E

# Saving data:
print "Saving Data"
name = 'rundat/Rot2D3C_FluxesAvg_'+str(date.today())+'.p'
pickle.dump(Data, open(str(name), 'wb'))
