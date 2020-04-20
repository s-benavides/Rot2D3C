import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import string
import os.path

# Plots energy as a function of time
# Assumes balance.txt is the output of an HD run
# Execute in ipython with '%run plot_energy.py'

# Scale for coloring:
def cscale(tau,taus):
    maxt = np.sqrt(taus[max(taus, key=taus.get)])
    mint = np.sqrt(taus[min(taus, key=taus.get)])
    return (np.sqrt(tau) - mint)/(maxt-mint)

# Making rule for reading nui,nun,mu
rule = string.maketrans('d', '0')

# Path to the data
num = input("how many runs to compare?")
runnames = []

# For saving:
svfig = raw_input("Save figures? (Y or N): ")
if svfig in ['Y']:
	svname = raw_input("Name for save: ")


for i in range(num):
        folder = raw_input("Folder name: ")
        runnames.append(folder)

omegas = dict([])
# For making colorbar
for jj,runname in enumerate(runnames):
  	path = '../'+runname+'/run/'
	omega = np.genfromtxt(path+'parameter.inp',comments='%',converters={0:  lambda val: float(val.translate(rule))},usecols=0)[19]
	omegas[runname]=omega

for jj,runname in enumerate(runnames):
  	path = '../'+runname+'/run/'

	params = np.genfromtxt(path+'parameter.inp',comments='%',converters={0:  lambda val: float(val.translate(rule))},usecols=0)

	rand= params[21]
	print('rand = ',rand)

	# Reads balance.txt
	t,en,enz,inj,injz,diss,dissz,hdiss,hdissz,coup,uf,ufz = np.loadtxt(path+'energy_bal.txt',unpack=True)

	if rand==3:
        	inj = inj/2. # RANDOM FORCING
        	injz = injz/2. # RANDOM FORCING

	# Averaging injection
	injtot = np.mean(inj+injz)

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
	if (len(runnames)==1):
		plt.figure(1)
		plt.title(runname)
		plt.plot(t,en,'--r',label='KE_2D')
		plt.plot(t,enz,'--g',label='KE_z')
		plt.plot(t,en+enz,'-k',label='KE_tot')
		
		plt.figure(2)
		plt.title(runname)
		plt.plot(t,hdiss,'--r',label='Hypo_2D')
		plt.plot(t,hdissz,'--g',label='Hypo_z')
		plt.plot(t,hdiss+hdissz,'-k',label='Hypo_tot')

		plt.figure(3)
		plt.title(runname)
		plt.plot(t,coup,'-b',label='Coupling')
                
                plt.figure(4)
                plt.title(runname)
                plt.plot(t,diss+dissz,'-k',label='disstot')
                plt.plot(t,diss,'--r',label='diss')
                plt.plot(t,dissz,'--g',label='dissz')

	else:
		plt.figure(1)
		plt.plot(t,en+enz,label=r"$\Omega = %.3f$" % omegas[runname],color = cm.copper(cscale(omegas[runname],omegas)))	
		#plt.plot(t,en+enz,label=runname)
	
		plt.figure(2)
		plt.plot(t,hdiss+hdissz,label=r"$\Omega = %.3f$" % omegas[runname],color = cm.copper(cscale(omegas[runname],omegas)))	
		#plt.plot(t,hdiss+hdissz,label = runname)
	
		plt.figure(3)
		plt.plot(t,coup,label=r"$\Omega = %.3f$" % omegas[runname],color = cm.copper(cscale(omegas[runname],omegas)))	
		#plt.plot(t,coup,label=runname)	

                plt.figure(4)
		plt.plot(t,diss+dissz,label=r"$\Omega = %.3f$" % omegas[runname],color = cm.copper(cscale(omegas[runname],omegas)))	
		#plt.plot(t,inj+injz,label=runname)	

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
plt.xlabel("Time")
plt.ylabel(r"$KE$")
plt.legend()
plt.tight_layout()
if svfig in ['Y']: plt.savefig('./figs/KE_'+svname+'.png')


plt.figure(2)
plt.xlabel("Time")
if (len(runnames)>1):
	plt.ylabel(r'Hypodissipation')
plt.legend()
plt.tight_layout()
if svfig in ['Y']: plt.savefig('./figs/Hypo_'+svname+'.png')

plt.figure(3)
if (len(runnames)>1):
#	plt.title(r"$|\phi|^2$")
	plt.ylabel(r"$\langle 2 \Omega \partial_y \psi v_z \rangle$")
plt.xlabel("Time")
plt.legend()
plt.tight_layout()
if svfig in ['Y']: plt.savefig('./figs/Coup_'+svname+'.png')

plt.figure(4)
if (len(runnames)>1):
	plt.ylabel(r"Viscous Dissipation")
plt.xlabel("Time")
plt.legend()
plt.tight_layout()
if svfig in ['Y']: plt.savefig('./figs/Diss_'+svname+'.png')


plt.show()
