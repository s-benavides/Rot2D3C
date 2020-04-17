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
    maxt = taus[max(taus, key=taus.get)]
    mint = taus[min(taus, key=taus.get)]
    return (tau - mint)/(maxt-mint)

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
	t,en,inj,den,hen = np.loadtxt(path+'u_bal.txt',unpack=True)
	t1,enphi,denphi,nlphi1,nlphi2 = np.loadtxt(path+'phi_bal.txt',unpack=True)
	t2,efk,e1,e2 = np.loadtxt(path+'misc.txt',unpack=True)
	t3,m1,m2,m3,m4 = np.loadtxt(path+'moments.txt',unpack=True)

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

	# Plots
	if (len(runnames)==1):
		plt.figure(1)
		plt.title(runname)
		plt.plot(t,en,'-k',label='KE')
		
		plt.figure(2)
		plt.title(runname)
		plt.plot(t,hen,'-k',label='Hypo')

		plt.figure(3)
		plt.title(runname)
		plt.plot(t1,enphi,'-k',label='E_phi')
		plt.plot(t1,nlphi2,'--r',label='omega*phi')
		plt.plot(t1,nlphi1,'-.b',label='phi^3')
		plt.plot(t1,denphi,'-g',label='Diss')
                
                plt.figure(4)
                plt.title(runname)
                plt.semilogy(t1,enphi,'-k',label='E_phi')
                plt.semilogy(t1,nlphi2,'--r',label='omega*phi')
                plt.semilogy(t1,nlphi1,'-.b',label='phi^3')
                plt.semilogy(t1,denphi,'-g',label='Diss')

                plt.figure(5)
                plt.title(runname)
                plt.plot(t3,m1,'-k',label='m1')
                plt.plot(t3,m2,'--r',label='m2')
                plt.plot(t3,m3,'-.b',label='m3')
                plt.plot(t3,m4,'-g',label='m4')

                plt.figure(6)
                plt.title(runname)
                plt.semilogy(t3,m1,'-k',label='m1')
                plt.semilogy(t3,m2,'--r',label='m2')
                plt.semilogy(t3,m3,'-.b',label='m3')
                plt.semilogy(t3,m4,'-g',label='m4')


	else:
		plt.figure(1)
		plt.plot(t,en,label=runname)
	
		plt.figure(2)
		plt.plot(t,hen,label = runname)
	
		plt.figure(3)
		plt.plot(t1,enphi,label=r"$\mu/\mu_c = %.3f$" % mus[runname],color = cm.copper(cscale(mus[runname],mus)))	

                plt.figure(4)
                plt.semilogy(t1,enphi,label=r"$\mu/\mu_c = %.3f$" % mus[runname],color = cm.copper(cscale(mus[runname],mus)))


	mufk = np.mean(efk)
	print('run: %s, mean ufk: %.3e' % (runname,mufk))
		
	print('run: %s, mean injtot: %f4' % (runname,injtot))
	
	Re_rms=np.sqrt(mufk/(nu*(kf)**(2*nn-1)))
	print('run: %s, Re_rms: %f4' % (runname,Re_rms))
	
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
	plt.ylabel(r"$\langle |\phi|^2 \rangle$")
plt.xlabel("Time")
plt.legend()
plt.tight_layout()
if svfig in ['Y']: plt.savefig('./figs/Enphi_'+svname+'.png')

plt.figure(4)
if (len(runnames)>1):
	plt.ylabel(r"$\langle |\phi|^2 \rangle$")
plt.xlabel("Time")
plt.legend()
plt.tight_layout()
if svfig in ['Y']: plt.savefig('./figs/Enphi_semilog_'+svname+'.png')


plt.show()
