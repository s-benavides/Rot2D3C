import numpy as np
import matplotlib.pyplot as plt

kf = 8#input('k_f = ')
rd = raw_input('Random forcing? Y or N: ')

t,en,inj,den,hen = np.loadtxt('./u_bal.txt',unpack=True)
tphi,enphi,denphi,nlphi,nlomega = np.loadtxt('./phi_bal.txt',unpack=True)
start=10#int(datain2[1]*10) # 10 comes from sstep/cstep
print('Start %s.'% (start))

###################################### Energy balance
diss = den[start:]+hen[start:]
inj = inj[start:]

if rd in ['Y']:
	inj = 0.5*inj #for random forcing
	print("Random")

dE = np.gradient(en[start:])/2.
dt = np.gradient(t[start:])
deriv = dE/dt

# Individual balance:
slope = inj - diss
diff = deriv- slope

dd = np.nanmean(diff)
ds = np.nanstd(diff,ddof=1)/float(np.sqrt(len(diff)))

print('Mean (dE/dt - Inj + Diss) = ',dd,'+/-',ds)

plt.figure(1)
plt.plot(diff,'-k',label=r"d E/ dt - Inj + Diss",lw=0.5)
#plt.plot(slope,'.b',label=r"Inj - Diss",alpha=0.5)
#plt.plot(deriv,'.r',label=r"dE/dt",alpha=0.5)
plt.plot(diff*0,'--r')
plt.legend()
plt.xlabel("Output (not time)")

################################## Phi balance
dE = np.gradient(enphi)/2.
dt = np.gradient(tphi)
deriv = dE/dt

# Individual balance:
slope = 0.01*enphi - 1.0*nlphi - denphi + nlomega
diff = deriv- slope

dd = np.nanmean(diff)
ds = np.nanstd(diff,ddof=1)/float(np.sqrt(len(diff)))

print('Mean (dPhi/dt - RHS) = ',dd,'+/-',ds)

plt.figure(2)
plt.title('Phi')
plt.plot(diff,'-k',label=r"d E/ dt - Inj + Diss",lw=0.5)
#plt.plot(slope,'.b',label=r"Inj - Diss",alpha=0.5)
#plt.plot(deriv,'.r',label=r"dE/dt",alpha=0.5)
plt.plot(diff*0,'--r')
plt.legend()
plt.xlabel("Output (not time)")


plt.show()
