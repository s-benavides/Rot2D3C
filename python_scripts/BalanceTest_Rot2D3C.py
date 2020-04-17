import numpy as np
import matplotlib.pyplot as plt

kf = 8#input('k_f = ')
rd = raw_input('Random forcing? Y or N: ')

path = './'

# Reads balance.txt
t,en,enz,inj,injz,diss,dissz,hdiss,hdissz,coup,uf,ufz = np.loadtxt(path+'energy_bal.txt',unpack=True)

start=10#int(datain2[1]*10) # 10 comes from sstep/cstep
print('Start %s.'% (start))

########################### Energy balance
### 2D energy ###
diss = diss[start:]+hdiss[start:]
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
plt.title("2D energy bal")
plt.legend()
plt.xlabel("Output (not time)")

### energy_z ###
t,en,enz,inj,injz,diss,dissz,hdiss,hdissz,coup,uf,ufz = np.loadtxt(path+'energy_bal.txt',unpack=True)
diss = dissz[start:]+hdissz[start:]
inj = injz[start:]

if rd in ['Y']:
        inj = 0.5*inj #for random forcing
        print("Random")

dE = np.gradient(enz[start:])/2.
dt = np.gradient(t[start:])
deriv = dE/dt

# Individual balance:
slope = inj - diss
diff = deriv- slope

dd = np.nanmean(diff)
ds = np.nanstd(diff,ddof=1)/float(np.sqrt(len(diff)))

print('Mean (dE/dt - Inj + Diss) = ',dd,'+/-',ds)

plt.figure(2)
plt.plot(diff,'-k',label=r"d E/ dt - Inj + Diss",lw=0.5)
#plt.plot(slope,'.b',label=r"Inj - Diss",alpha=0.5)
#plt.plot(deriv,'.r',label=r"dE/dt",alpha=0.5)
plt.plot(diff*0,'--r')
plt.title("energy_z bal")
plt.legend()
plt.xlabel("Output (not time)")

### total energy (if omega=/=0) ###
t,en,enz,inj,injz,diss,dissz,hdiss,hdissz,coup,uf,ufz = np.loadtxt(path+'energy_bal.txt',unpack=True)
diss = diss[start:]+hdiss[start:]+dissz[start:]+hdissz[start:]
inj = inj[start:]+injz[start:]

if rd in ['Y']:
        inj = 0.5*inj #for random forcing
        print("Random")

dE = np.gradient(en[start:]+enz[start:])/2.
dt = np.gradient(t[start:])
deriv = dE/dt

# Individual balance:
slope = inj - diss
diff = deriv- slope

dd = np.nanmean(diff)
ds = np.nanstd(diff,ddof=1)/float(np.sqrt(len(diff)))

print('Mean (dE/dt - Inj + Diss) = ',dd,'+/-',ds)

plt.figure(3)
plt.plot(diff,'-k',label=r"d E/ dt - Inj + Diss",lw=0.5)
#plt.plot(slope,'.b',label=r"Inj - Diss",alpha=0.5)
#plt.plot(deriv,'.r',label=r"dE/dt",alpha=0.5)
plt.plot(diff*0,'--r')
plt.title("total energy bal")
plt.legend()
plt.xlabel("Output (not time)")

################################## Enstrophy balance
t,en,inj,diss,hdiss = np.loadtxt(path+'enstrophy_bal.txt',unpack=True)

diss = diss[start:]+hdiss[start:]
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

print('Mean (dZ/dt - Inj + Diss) = ',dd,'+/-',ds)

plt.figure(4)
plt.plot(diff,'-k',label=r"d Z/ dt - Inj + Diss",lw=0.5)
#plt.plot(slope,'.b',label=r"Inj - Diss",alpha=0.5)
#plt.plot(deriv,'.r',label=r"dE/dt",alpha=0.5)
plt.plot(diff*0,'--r')
plt.title("enstrophy bal")
plt.legend()
plt.xlabel("Output (not time)")

################################## Helicity balance
t,en,inj,injz,diss,dissz,hdiss,hdissz = np.loadtxt(path+'helicity_bal.txt',unpack=True)

diss = diss[start:]+hdiss[start:]+dissz[start:]+hdissz[start:]
inj = inj[start:]+injz[start:]

# Temporary solution to sign error
inj = -inj

if rd in ['Y']:
        inj = 0.5*inj #for random forcing
        print("Random")

dE = np.gradient(en[start:])
dt = np.gradient(t[start:])
deriv = dE/dt

# Individual balance:
slope = inj - diss
diff = deriv- slope

dd = np.nanmean(diff)
ds = np.nanstd(diff,ddof=1)/float(np.sqrt(len(diff)))

print('Mean (dH/dt - Inj + Diss) = ',dd,'+/-',ds)

plt.figure(5)
plt.plot(diff,'-k',label=r"d H/ dt - Inj + Diss",lw=0.5)
#plt.plot(slope,'.b',label=r"Inj - Diss",alpha=0.5)
#plt.plot(deriv,'.r',label=r"dE/dt",alpha=0.5)
plt.plot(diff*0,'--r')
plt.title("Helicity bal")
plt.legend()
plt.xlabel("Output (not time)")

plt.show()
