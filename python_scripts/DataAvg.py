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
avglist = sorted(glob.glob('rundat/AvgTimemu*.txt'))

runnames = []
for file in avglist:
  run = file.split('rundat/AvgTime')[1]
  run = run.split('.txt')[0]
  runnames.append(run)

Data = dict([])
for i,run in enumerate(runnames):
	Data_E = dict([])
	print("Working on run %s " % run)
	path = '../'+run+'/run/'

	params = np.genfromtxt(path+'parameter.inp',comments='%',converters={0:  lambda val: float(val.translate(rule))},usecols=0)

        rand= params[18]

	sstep = params[3]

	cstep = params[2]


        # Reads balance.txt
        t,en,inj,den,hen = np.loadtxt(path+'u_bal.txt',unpack=True)
        t1,enphi,denphi,nlphi1,nlphi2 = np.loadtxt(path+'phi_bal.txt',unpack=True)
        t2,efk,e1,e2 = np.loadtxt(path+'misc.txt',unpack=True)

        if rand==3:
                inj = inj/2. # RANDOM FORCING

        # Averaging injection
        injtot = np.mean(inj)

        # nu,kf
        nu = float(params[10])
        hnu = float(params[11])
        nn = float(params[12])
        mm = float(params[13])
        kdn = float(params[8])
        kup = float(params[9])
        kf = (kdn+kup)/2.

        # mu
        mu = float(params[14])
        cphi = float(params[15])
        dphi = float(params[16])
        fphi = float(params[22])

	# Average start indices
	[start,start_fl] = np.loadtxt('rundat/AvgTime'+run+'.txt')

	start = int(start)
	start_fl = int(start)

	# list of observables to average:
	olist = {'en':en,'inj':inj,'den':den,'hen':hen,'enphi':enphi,'denphi':denphi,'nlphi1':nlphi1,'nlphi2':nlphi2,'efk':efk,'e1':e1,'e2':e2}
	
        # AVERAGING
	Data_E = dict([])
	for obs in olist:
		avg = np.nanmean(olist[obs][start:])
		err = np.std(olist[obs][start:])/np.sqrt(len(olist[obs][start:]))
		Data_E[obs]=[avg,err]

        # Some calculations
#        Re_rms=np.sqrt(avg_ufk)/(nu*(kf)**(2*hek-1))
        Re_rms=np.sqrt(Data_E['efk'][0])/(nu*(kf)**(2*nn-1))


	Data_E['Re_rms'] = Re_rms
	Data_E['nu'] = nu
	Data_E['hnu'] = hnu
	Data_E['mu'] = mu
	Data_E['nn'] = nn
	Data_E['mm'] = mm
	Data_E['kdn'] = kdn
	Data_E['kup'] = kup
	Data_E['kf'] = kf
	Data_E['cphi'] = cphi
	Data_E['dphi'] = dphi
	Data_E['fphi'] = fphi
		
	Data[run] = Data_E
	
# Saving data:
print "Saving Data"
name = 'rundat/DataAvg_'+str(date.today())+'.p'
pickle.dump(Data, open(str(name), 'wb'))	
