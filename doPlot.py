import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp
from matplotlib.ticker import MaxNLocator

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


nTicks = 5

f = open('posterior.sizes', 'r')
dat = f.read()
dat = dat.strip()
[npar, nstep] = map(int, dat.split())

ch = np.fromfile('test.extract', dtype=np.float64).reshape((npar,nstep), order="FORTRAN")

length = len(ch[0,:])
#ch = ch[:,length/2:]
ch = ch[:,-5000:]

plt.close('all')

fig1 = plt.figure(num=1, figsize=(17,10))
plt.clf()

loop = 1
labels = [r'$\alpha_B$',r'$\beta_B$',r'$\alpha_\mu$',r'$\beta_\mu$',r'$\alpha_f$',r'$\beta_f$']
for i in range(4):
	ax = fig1.add_subplot(2,5,loop)
	ax.plot(ch[i,:])
	ax.set_xlabel('Iteration')
	ax.set_ylabel(labels[i])
	ax.xaxis.set_major_locator(MaxNLocator(nTicks))
	loop += 1
	ax = fig1.add_subplot(2,5,loop)
	ax.hist(ch[i,:])
	ax.set_xlabel(labels[i])
	ax.set_ylabel('p('+labels[i]+'|D)')
	ax.xaxis.set_major_locator(MaxNLocator(nTicks))
	loop += 1
	if ((i+1) % 2 == 0):
		loop += 1
	
# Magnetic field strength
B = np.linspace(0.1,500,100)
pB = np.zeros(100)
alpha = ch[0,:]
beta = ch[1,:]
pB = IGAvgPrior(B, alpha, beta)
ax = fig1.add_subplot(2,5,5)
ax.plot(B,pB)
ax.set_xlabel(r'B [G]')
ax.set_ylabel(r'p(B)')
ax.xaxis.set_major_locator(MaxNLocator(nTicks))

# Inclination
left = -1.0
right = 1.0
mu = np.linspace(left + 1e-2,right - 1e-2,100)
pmu = np.zeros(100)
alpha = ch[2,:]
beta = ch[3,:]
pmu = betaAvgPrior(mu, alpha, beta, left, right)
ax = fig1.add_subplot(2,5,10)
ax.plot(mu,pmu)
ax.set_xlabel(r'$\mu$')
ax.set_ylabel(r'p($\mu$)')
ax.xaxis.set_major_locator(MaxNLocator(nTicks))

# Filling factor
#left = 0.0
#right = 1.0
#f = np.linspace(left + 1e-4, right - 1e-4, 100)
#pf = np.zeros(100)
#alpha = ch[4,:]
#beta = ch[5,:]
#pf = betaAvgPrior(f, alpha, beta, left, right)
#ax = fig1.add_subplot(3,5,15)
#ax.plot(f,pf)
#ax.set_xlabel('f')
#ax.set_ylabel('p(f)')
#ax.xaxis.set_major_locator(MaxNLocator(6))

fig1.tight_layout()
