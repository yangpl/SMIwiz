import matplotlib.pyplot as plt
import numpy as np

#-----------------------------------------------
z, x, y, azimuth, dip, ireceiver, isrc = np.loadtxt('sources.txt', skiprows=0, unpack=True)

batch1 = np.loadtxt('batch1', skiprows=0, unpack=True)
batch2 = np.loadtxt('batch2', skiprows=0, unpack=True)
batch3 = np.loadtxt('batch3', skiprows=0, unpack=True)
batch4 = np.loadtxt('batch4', skiprows=0, unpack=True)
batch5 = np.loadtxt('batch5', skiprows=0, unpack=True)
batch6 = np.loadtxt('batch6', skiprows=0, unpack=True)
batch7 = np.loadtxt('batch7', skiprows=0, unpack=True)
batch8 = np.loadtxt('batch8', skiprows=0, unpack=True)


plt.figure(figsize=(7, 7))
plt.scatter(x=x, y=y, s=50, c=z, cmap='rainbow')

for ishot in batch1:
    idx = (isrc==ishot)
    plt.scatter(x=x[idx], y=y[idx], s=50, color='red')
plt.scatter(x=x[idx], y=y[idx], s=50, color='red', label='batch1')
for ishot in batch2:
    idx = (isrc==ishot)
    plt.scatter(x=x[idx], y=y[idx], s=50, color='green')
plt.scatter(x=x[idx], y=y[idx], s=50, color='green', label='batch2')
for ishot in batch3:
    idx = (isrc==ishot)
    plt.scatter(x=x[idx], y=y[idx], s=50, color='blue')
plt.scatter(x=x[idx], y=y[idx], s=50, color='blue', label='batch3')
for ishot in batch4:
    idx = (isrc==ishot)
    plt.scatter(x=x[idx], y=y[idx], s=50, color='y')
plt.scatter(x=x[idx], y=y[idx], s=50, color='y', label='batch4')
for ishot in batch5:
    idx = (isrc==ishot)
    plt.scatter(x=x[idx], y=y[idx], s=50, color='k')
plt.scatter(x=x[idx], y=y[idx], s=50, color='k', label='batch5')
for ishot in batch6:
    idx = (isrc==ishot)
    plt.scatter(x=x[idx], y=y[idx], s=50, color='m')
plt.scatter(x=x[idx], y=y[idx], s=50, color='m', label='batch6')
for ishot in batch7:
    idx = (isrc==ishot)
    plt.scatter(x=x[idx], y=y[idx], s=50, color='grey')
plt.scatter(x=x[idx], y=y[idx], s=50, color='grey', label='batch7')
for ishot in batch8:
    idx = (isrc==ishot)
    plt.scatter(x=x[idx], y=y[idx], s=50, color='purple')
plt.scatter(x=x[idx], y=y[idx], s=50, color='purple', label='batch8')

plt.legend(bbox_to_anchor=(0.93, 0.5))    
plt.xlabel('X(m)')
plt.ylabel('Y(m)')


plt.savefig('scatter.png', bbox_inches='tight')
plt.show()

