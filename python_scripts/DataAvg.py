import numpy as np
import matplotlib.pyplot as plt
#from matplotlib import rc
import pickle
import glob as glob
import string
from datetime import date
import bunch_err
import field_calc

# Making rule for reading nui,nun,mu
rule = string.maketrans('d', '0')

# Get all the AvgTimeO*B*.txt files
avglist = sorted(glob.glob('rundat/AvgTimeO*.txt'))

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
	
	rand = params[21]
	
	# Reads balance.txt
	t,en,enz,inj,injz,diss,dissz,hdiss,hdissz,coup,uf,ufz = np.loadtxt(path+'energy_bal.txt',unpack=True)
	t,enst,injenst,dissenst,hdissenst = np.loadtxt(path+'enstrophy_bal.txt',unpack=True)	
	t,hel,injhel1,injhel2,disshel1,disshel2,hdisshel1,hdisshel2 = np.loadtxt(path+'helicity_bal.txt',unpack=True)
	
	if rand==3:
	        inj = inj/2. # RANDOM FORCING
	        injz = injz/2. # RANDOM FORCING
	        injenst = injenst/2. # RANDOM FORCING
	        injhel1 = injhel1/2. # RANDOM FORCING
	        injhel2 = injhel2/2. # RANDOM FORCING
	
	# nu,kf
	nu = float(params[11])
	hnu = float(params[12])
	nn = float(params[13])
	mm = float(params[14])
	kdn = float(params[9])
	kup = float(params[10])
	kf = (kdn+kup)/2.
	
	# nuv, etc
	nuv = float(params[15])
	hnuv = float(params[16])
	nnv = float(params[17])
	mmv = float(params[18])
	omega = float(params[19])

	# Average start indices
	[start,start_fl,err_ind] = np.loadtxt('rundat/AvgTime'+run+'.txt')

	start = int(start)
	start_fl = int(start)
	err_ind = int(err_ind)

	# list of observables to average:
	olist = {'en':en,'inj':inj,'diss':diss,'hdiss':hdiss,'enz':enz,'injz':injz,'dissz':dissz,'hdissz':hdissz,'coup':coup,'uf':uf,'ufz':ufz,
		'enst':enst,'injenst':injenst,'dissenst':dissenst,'hdissenst':hdissenst,
		'hel':hel,'injhel1':injhel1,'disshel1':disshel1,'hdisshel1':hdisshel1,'injhel2':injhel2,'disshel2':disshel2,'hdisshel2':hdisshel2,
		}
	
        # AVERAGING
	Data_E = dict([])
	for obs in olist:
		avg = np.nanmean(olist[obs][start:])
		err = bunch_err.bunch_err(olist[obs][start:],err_ind=err_ind)
		Data_E[obs]=[avg,err]

        # Some calculations
#        Re_rms=np.sqrt(avg_ufk)/(nu*(kf)**(2*hek-1))
        Re_rms=np.sqrt(Data_E['uf'][0]+Data_E['ufz'][0])/(nu*(kf)**(2*nn-1))
        Re_2D_rms=np.sqrt(Data_E['uf'][0])/(nu*(kf)**(2*nn-1))
        Re_z_rms=np.sqrt(Data_E['ufz'][0])/(nuv*(kf)**(2*nnv-1))

	# Components of velocity
        reso = 512
	tf = np.loadtxt('../'+run+'/run/time_field.txt')
        outnum = str(int(tf[-1][0])) #raw_input("out num? ") #sys.argv[1]
        outnum ="{:0>3s}".format(outnum)	

	otypes = ['vx','vy']
	for otype in otypes:
        	out = field_calc.field_calc(run,otype,outnum,reso=reso,num_files=16)
                mean = 0.0
                mean = np.sum(np.abs(out)**2)/float(reso)**2
                Data_E[otype]=mean

	Data_E['rand'] = rand
	Data_E['Re_rms'] = Re_rms
	Data_E['Re_2D_rms'] = Re_2D_rms
	Data_E['Re_z_rms'] = Re_z_rms
	Data_E['nu'] = nu
	Data_E['hnu'] = hnu
	Data_E['nn'] = nn
	Data_E['mm'] = mm
	Data_E['kdn'] = kdn
	Data_E['kup'] = kup
	Data_E['kf'] = kf
        Data_E['nuv'] = nuv
        Data_E['hnuv'] = hnuv
        Data_E['nnv'] = nnv
        Data_E['mmv'] = mmv
        Data_E['omega'] = omega
	
	Data[run] = Data_E
	
# Saving data:
print "Saving Data"
name = 'rundat/Rot2D3C_DataAvg_'+str(date.today())+'.p'
pickle.dump(Data, open(str(name), 'wb'))	
