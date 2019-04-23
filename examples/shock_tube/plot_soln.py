import matplotlib.pyplot as plt
import numpy as np
        
data_file = 'x.dat'
data      = np.loadtxt(data_file, unpack=False)

plt.figure(1)
plt.xlabel('x',fontsize=20)
plt.ylabel('Density',fontsize=20) 
plt.plot(data[:,0],data[:,1],lw=2)
plt.axis('tight')
plt.draw()
plt.savefig('density.eps', format='eps')
plt.clf()

plt.figure(1)
plt.xlabel('x',fontsize=20)
plt.ylabel('Velocity',fontsize=20) 
plt.plot(data[:,0],data[:,2],lw=2)
plt.axis('tight')
plt.draw()
plt.savefig('vel_x.eps', format='eps')
plt.clf()

plt.figure(1)
plt.xlabel('x',fontsize=20)
plt.ylabel('Pressure',fontsize=20) 
plt.plot(data[:,0],data[:,3],lw=2)
plt.axis('tight')
plt.draw()
plt.savefig('pressure.eps', format='eps')
plt.clf()


