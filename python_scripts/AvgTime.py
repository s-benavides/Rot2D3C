import numpy as np
import matplotlib.pyplot as plt
import string
import bunch_err
# Plots energy as a function of time
# Assumes balance.txt is the output of an HD run
# Execute in ipython with '%run plot_energy.py'

# Making rule for reading nui,nun,mu
rule = string.maketrans('d', '0')

# Path to the data
runname = raw_input("What folder name? ")

path = '../'+runname+'/run/'

params = np.genfromtxt(path+'parameter.inp',comments='%',converters={0:  lambda val: float(val.translate(rule))},usecols=0)

rand = params[21]

sstep = params[3]

cstep = params[2]

# Reads balance.txt
t,en,enz,inj,injz,diss,dissz,hdiss,hdissz,coup,uf,ufz = np.loadtxt(path+'energy_bal.txt',unpack=True)

if rand==3:
	inj = inj/2. # RANDOM FORCING
	injz = injz/2. # RANDOM FORCING

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

# Plots
plt.figure(1)
plt.title(runname)
plt.plot(en,'--r',label=r'$E_{2D}$')
plt.plot(enz,'--g',label=r'$E_{z}$')
plt.plot(en+enz,'-k',label=r'$E_{tot}$')
	
plt.figure(2)
plt.title(runname)
plt.plot(hdiss,'--r',label = 'Hypo_2D')
plt.plot(hdissz,'--g',label = 'Hypo_z')
plt.plot(hdiss+hdissz,'-k',label = 'Hypo')

mufk = np.mean(uf)
print('run: %s, mean ufk: %.3e' % (runname,mufk))		
print('run: %s, mean inj: %f4' % (runname,np.mean(inj)))
Re_rms=np.sqrt(mufk)/(nu*(kf)**(2*nn-1))
print('run: %s, Re_rms: %f4' % (runname,Re_rms))
	
mufk = np.mean(ufz)
print('run: %s, mean ufz: %.3e' % (runname,mufk))
print('run: %s, mean inz: %f4' % (runname,np.mean(injz)))
Re_rms=np.sqrt(mufk)/(nuv*(kf)**(2*nnv-1))
print('run: %s, Re_rms_z: %f4' % (runname,Re_rms))


plt.figure(1)
plt.xlabel("Output #")
plt.legend()
plt.tight_layout()

plt.figure(2)
plt.xlabel("Output #")
#plt.ylim(0,1.0)
plt.legend()
plt.tight_layout()

plt.show()

start = input('Enter x-axis value you want to start averaging: ')
start_fl = start/float(sstep/cstep)        #starting flux and spectra number
print('Error for hypodiss')
bunch_err.bunch_err(hdiss[start:])
err_ind=input('Enter iteration number for error calc: ')

np.savetxt('rundat/AvgTime'+runname+'.txt',[start,start_fl,err_ind], delimiter ='   ')

