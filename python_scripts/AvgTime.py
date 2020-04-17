import numpy as np
import matplotlib.pyplot as plt
import string

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

# Averaging injection
injtot = np.mean(inj) + np.mean(injz)

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

# Plots
plt.figure(1)
plt.title(runname)
plt.plot(en,'-k',label='E')
	
plt.figure(2)
plt.title(runname)
plt.plot(hen,'-k',label = 'Hypo')

plt.figure(3)
plt.title(runname)
plt.plot(enphi,'-k',label = 'enphi')

plt.figure(4)
plt.title(runname)
plt.semilogy(enphi,'-k',label = 'enphi')

mufk = np.mean(efk)
print('run: %s, mean ufk: %.3e' % (runname,mufk))		
print('run: %s, mean inj: %f4' % (runname,injtot))
Re_rms=np.sqrt(mufk)/(nu*(kf)**(2*nn-1))
print('run: %s, Re_rms: %f4' % (runname,Re_rms))
	
plt.figure(1)
plt.xlabel("Output #")
plt.ylabel(r"$E_{tot}$")
plt.legend()
plt.tight_layout()

plt.figure(2)
plt.xlabel("Output #")
plt.ylabel("Hypovisc")
#plt.ylim(0,1.0)
plt.legend()
plt.tight_layout()

plt.figure(3)
plt.xlabel("Output #")
plt.ylabel(r"$|\phi|^2$")
#plt.ylim(0,1.0)
plt.legend()
plt.tight_layout()

plt.figure(4)
plt.xlabel("Output #")
plt.ylabel(r"$|\phi|^2$")
#plt.ylim(0,1.0)
plt.legend()
plt.tight_layout()

plt.show()

start = input('Enter x-axis value you want to start averaging: ')
start_fl = start/float(sstep/cstep)        #starting flux and spectra number

np.savetxt('rundat/AvgTime'+runname+'.txt',[start,start_fl], delimiter ='   ')

