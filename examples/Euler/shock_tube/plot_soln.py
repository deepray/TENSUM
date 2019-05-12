#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import sys
        
data_file = sys.argv[1]
ref_file  = 'ref_soln.dat'
data      = np.loadtxt(data_file, unpack=False)
ref_data  = np.loadtxt(ref_file, unpack=False)

plt.figure(1)
plt.xlabel('x',fontsize=20)
plt.ylabel('Density',fontsize=20) 
plt.plot(ref_data[:,0],ref_data[:,1],lw=2,label='reference',color='black')
plt.plot(data[:,0],data[:,1],lw=2,label='numerical',color='red')
plt.legend(loc='upper right',fontsize=15)
plt.xlim(0,1)
plt.ylim(0,1.2)
plt.draw()
plt.savefig('density.pdf',dvi=300,bbox_inches='tight')
plt.clf()

plt.figure(1)
plt.xlabel('x',fontsize=20)
plt.ylabel('Velocity',fontsize=20) 
plt.plot(ref_data[:,0],ref_data[:,2],lw=2,label='reference',color='black')
plt.plot(data[:,0],data[:,2],lw=2,label='numerical',color='red')
plt.legend(loc='lower center',fontsize=15)
plt.xlim(0,1)
plt.ylim(0,1.5)
plt.draw()
plt.savefig('vel_x.pdf',dvi=300,bbox_inches='tight')
plt.clf()

plt.figure(1)
plt.xlabel('x',fontsize=20)
plt.ylabel('Pressure',fontsize=20) 
plt.plot(ref_data[:,0],ref_data[:,3],lw=2,label='reference',color='black')
plt.plot(data[:,0],data[:,3],lw=2,label='numerical',color='red')
plt.legend(loc='upper right',fontsize=15)
plt.xlim(0,1)
plt.ylim(0,1.2)
plt.draw()
plt.savefig('pressure.pdf',dvi=300,bbox_inches='tight')
plt.clf()


