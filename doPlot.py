import numpy as np
import matplotlib.pyplot as plt

f = open('posterior.sizes', 'r')
dat = f.read()
dat = dat.strip()
[npar, nstep] = map(int, dat.split())

ch = np.fromfile('posterior', dtype=np.float64).reshape((npar,nstep), order="FORTRAN")

fig1 = plt.figure()

ax = fig1.add_subplot(111)
ax.plot(ch[0,:] * np.cos(np.pi*ch[1,:]/180.))