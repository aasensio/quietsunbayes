import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp

def betaAvgPrior(x, alpha, beta, left, right):
	Beta = sp.beta(alpha, beta)
	pf = np.zeros(len(x))
	for i in range(len(x)):		
		y = (right-left)**(1.0-alpha-beta) / Beta * (x[i] - left)**(alpha-1.0) * (right - x[i])**(beta-1.0)
		pf[i] = np.mean(y)
	return pf

def IGAvgPrior(x, alpha, beta):
	pf = np.zeros(len(x))
	for i in range(len(x)):		
		y = beta**alpha / sp.gamma(alpha) * x[i]**(-(alpha+1.0)) * np.exp(-beta/x[i])
		pf[i] = np.mean(y)
	return pf
	
f = open('posterior.sizes', 'r')
dat = f.read()
dat = dat.strip()
[npar, nstep] = map(int, dat.split())

ch = np.fromfile('posterior', dtype=np.float64).reshape((npar,nstep), order="FORTRAN")

length = len(ch[0,:])
ch = ch[:,length/2:]

fig1 = plt.figure(1)
plt.clf()

for i in range(8):
	ax = fig1.add_subplot(4,2,i+1)
	ax.plot(ch[i,:])
	
fig2 = plt.figure(2)
plt.clf()

# Magnetic field strength
B = np.linspace(0.1,2500,100)
pB = np.zeros(100)
alpha = ch[0,:]
beta = ch[1,:]
pB = IGAvgPrior(B, alpha, beta)
ax = fig2.add_subplot(2,2,1)
ax.plot(B,pB)

# Inclination
left = -1.0
right = 1.0
mu = np.linspace(left,right,100)
pmu = np.zeros(100)
alpha = ch[2,:]
beta = ch[3,:]
pmu = betaAvgPrior(mu, alpha, beta, left, right)
ax = fig2.add_subplot(2,2,2)
ax.plot(mu,pmu)

# Filling factor
left = 0.0
right = 1.0
f = np.linspace(left, right, 100)
pf = np.zeros(100)
alpha = ch[4,:]
beta = ch[5,:]
pf = betaAvgPrior(f, alpha, beta, left, right)
ax = fig2.add_subplot(2,2,3)
ax.plot(f,pf)

# Azimuth
left = 0.0
right = np.pi
phi = np.linspace(left, right, 100)
pphi = np.zeros(100)
alpha = ch[6,:]
beta = ch[7,:]
pphi = betaAvgPrior(phi, alpha, beta, left, right)
ax = fig2.add_subplot(2,2,4)
ax.plot(phi,pphi)