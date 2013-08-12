import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp

def betaAvgPrior(x, alpha, beta, left, right):
	Beta = sp.beta(alpha, beta)
	pf = np.zeros(len(x))
	for i in range(len(x)):
		ylog = ( (1.0-alpha-beta) * np.log(right-left) - (sp.gammaln(alpha) + sp.gammaln(beta) - sp.gammaln(alpha+beta)) +
			(alpha-1.0) * np.log(x[i] - left) + (beta-1.0) * np.log(right - x[i]) )		
		pf[i] = np.mean(np.exp(ylog))
	return pf

def IGAvgPrior(x, alpha, beta):
	pf = np.zeros(len(x))
	for i in range(len(x)):		
		logy = alpha * np.log(beta) - sp.gammaln(alpha) + (-(alpha+1.0)) * np.log(x[i]) -beta/x[i]
		pf[i] = np.mean(np.exp(logy))
	return pf
	
f = open('posterior.sizes', 'r')
dat = f.read()
dat = dat.strip()
[npar, nstep] = map(int, dat.split())

ch = np.fromfile('test.extract', dtype=np.float64).reshape((npar,nstep), order="FORTRAN")

length = len(ch[0,:])
ch = ch[:,length/2:]

fig1 = plt.figure(1, figsize=(15,10))
plt.clf()

loop = 1
for i in range(6):
	ax = fig1.add_subplot(3,5,loop)
	ax.plot(ch[i,:])	
	loop += 1
	ax = fig1.add_subplot(3,5,loop)
	ax.hist(ch[i,:])
	loop += 1
	if ((i+1) % 2 == 0):
		loop += 1
	
# Magnetic field strength
B = np.linspace(0.1,1200,100)
pB = np.zeros(100)
alpha = ch[0,:]
beta = ch[1,:]
pB = IGAvgPrior(B, alpha, beta)
ax = fig1.add_subplot(3,5,5)
ax.plot(B,pB)

# Inclination
left = -1.0
right = 1.0
mu = np.linspace(left + 1e-4,right - 1e-4,100)
pmu = np.zeros(100)
alpha = ch[2,:]
beta = ch[3,:]
pmu = betaAvgPrior(mu, alpha, beta, left, right)
ax = fig1.add_subplot(3,5,10)
ax.plot(mu,pmu)

# Filling factor
left = 0.0
right = 1.0
f = np.linspace(left + 1e-4, right - 1e-4, 100)
pf = np.zeros(100)
alpha = ch[4,:]
beta = ch[5,:]
pf = betaAvgPrior(f, alpha, beta, left, right)
ax = fig1.add_subplot(3,5,15)
ax.plot(f,pf)

fig1.tight_layout()
