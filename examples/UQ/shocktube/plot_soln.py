#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import sys

Sample_start = int(sys.argv[1])
Sample_end   = int(sys.argv[2])
tstamp       = int(sys.argv[3])
 
# Loading data        
for s in range(Sample_start,Sample_end+1):   
    data_file = 'SAMPLE_'+str(s)+'/line_data_'+str(tstamp)+'.dat'     
    data      = np.loadtxt(data_file, unpack=False)
    
    if(s==Sample_start):
        x     = data[:,0:1];
        rho   = data[:,1:2];
        vel   = data[:,2:3];
        pre   = data[:,3:4];
    else:
        rho   = np.hstack((rho,data[:,1:2]));
        vel   = np.hstack((vel,data[:,2:3]));
        pre   = np.hstack((pre,data[:,3:4]));  
               
# Plotting all samples
for s in range(Sample_start,Sample_end+1): 
    plt.figure(1)
    plt.xlabel('x',fontsize=20)
    plt.ylabel('Density',fontsize=20) 
    #plt.plot(ref_data[:,0],ref_data[:,1],lw=2,label='reference',color='black')
    plt.plot(x,rho[:,s],lw=2,label='numerical')
    #plt.legend(loc='upper right',fontsize=15)
    plt.xlim(0,1)
    plt.ylim(0,1.2)
    plt.draw()
    

    plt.figure(2)
    plt.xlabel('x',fontsize=20)
    plt.ylabel('Velocity',fontsize=20) 
    #plt.plot(ref_data[:,0],ref_data[:,2],lw=2,label='reference',color='black')
    plt.plot(x,vel[:,s],lw=2,label='numerical')
    #plt.legend(loc='lower center',fontsize=15)
    plt.xlim(0,1)
    plt.ylim(-.1,1.6)
    plt.draw()


    plt.figure(3)
    plt.xlabel('x',fontsize=20)
    plt.ylabel('Pressure',fontsize=20) 
    #plt.plot(ref_data[:,0],ref_data[:,3],lw=2,label='reference',color='black')
    plt.plot(x,pre[:,s],lw=2,label='numerical')
    #plt.legend(loc='upper right',fontsize=15)
    plt.xlim(0,1)
    plt.ylim(0,1.2)
    plt.draw()

plt.figure(1)
plt.savefig('density_samples.pdf',dvi=300,bbox_inches='tight')

plt.figure(2)
plt.savefig('vel_x_samples.pdf',dvi=300,bbox_inches='tight')

plt.figure(3)
plt.savefig('pressure_samples.pdf',dvi=300,bbox_inches='tight')


# Evaluating pointwise statistics
# Density
pt_mean = np.mean(rho,axis=1);
pt_std  = np.std(rho,axis=1);
plt.figure(4)
plt.xlabel('x',fontsize=20)
plt.ylabel('Density',fontsize=20)
plt.plot(x,pt_mean,'-',lw=2,label='mean',color='blue')
plt.plot(x,pt_mean+pt_std,'--',lw=2,label='mean+std',color='red')
plt.plot(x,pt_mean-pt_std,'--',lw=2,label='mean-std',color='green')
plt.xlim(0,1)
plt.ylim(0,1.2)
plt.legend(loc='upper right',fontsize=15)
plt.draw()
plt.savefig('density_stats.pdf',dvi=300,bbox_inches='tight')
plt.clf()

# Velocity
pt_mean = np.mean(vel,axis=1);
pt_std  = np.std(vel,axis=1);
plt.figure(4)
plt.xlabel('x',fontsize=20)
plt.ylabel('Velocity',fontsize=20)
plt.plot(x,pt_mean,'-',lw=2,label='mean',color='blue')
plt.plot(x,pt_mean+pt_std,'--',lw=2,label='mean+std',color='red')
plt.plot(x,pt_mean-pt_std,'--',lw=2,label='mean-std',color='green')
plt.legend(loc='lower center',fontsize=15)
plt.xlim(0,1)
plt.ylim(-0.1,1.6)
plt.draw()
plt.savefig('vel_x_stats.pdf',dvi=300,bbox_inches='tight')
plt.clf()

# Pressure
pt_mean = np.mean(pre,axis=1);
pt_std  = np.std(pre,axis=1);
plt.figure(4)
plt.xlabel('x',fontsize=20)
plt.ylabel('Pressure',fontsize=20)
plt.plot(x,pt_mean,'-',lw=2,label='mean',color='blue')
plt.plot(x,pt_mean+pt_std,'--',lw=2,label='mean+std',color='red')
plt.plot(x,pt_mean-pt_std,'--',lw=2,label='mean-std',color='green')
plt.legend(loc='upper right',fontsize=15)
plt.xlim(0,1)
plt.ylim(0,1.2)
plt.draw()
plt.savefig('pressure_stats.pdf',dvi=300,bbox_inches='tight')
plt.clf()





