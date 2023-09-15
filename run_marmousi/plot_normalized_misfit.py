import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['axes.xmargin'] = 0

import numpy as np
import cmath
import os

#-----------------------------------------------
emf = np.loadtxt('iterate.txt', skiprows=8, unpack=True)
it = emf[0,:]
fk = emf[1,:]
fk_f0 = emf[2,:]
gk = emf[3,:]



plt.plot(it, fk_f0, 'k')
plt.xlabel('# iteration k')
plt.ylabel('$J_k/J_0$')
#plt.title('Normalized misfit')

plt.xticks(np.arange(0, 50, 5))
plt.savefig('normalized_misfit.png', bbox_inches='tight')
plt.show()


