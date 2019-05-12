#!/usr/bin/env python

import sys, glob
from math import sqrt, pi, atan, cos,sin
import numpy as np
import os
import matplotlib.pyplot as plt

Base_folder = sys.argv[1]
npart       = int(sys.argv[2])
surf_tag    = [int(x) for x in sys.argv[3:]]

assert(len(surf_tag) == 2)

# Farfield parameters from param.in
T     = 1.0
p     = 1.0
M     = 0.729
gam   = 1.4
gasC  = 1.0
aoa   = 2.31/180.0*pi 


# Forces
force_file  = Base_folder + '/force.dat'
fdata       = np.loadtxt(force_file, unpack=False)
f_conv      = fdata[-1,2::] # The final converged forces

ca    = cos(aoa)
sa    = sin(aoa)
uinf  = M*sqrt(gam*gasC*T)*ca
vinf  = M*sqrt(gam*gasC*T)*sa
KE    = 0.5*sqrt(uinf*uinf + vinf*vinf)*(gasC*T/p)

# Lift and drag coefficients due to pressure
cdp = (f_conv[0]*ca + f_conv[1]*sa)/KE
clp = (-f_conv[0]*sa + f_conv[1]*ca)/KE
print("Lift (pressure) = %.6f, Drag (pressure) = %.6f" %(clp,cdp))


# Cp data
# NOTE: By construction, each boundary face belongs to only one partition
#       thus there should be no duplicates. However they may not be in the 
#       right in order.    



for tag in surf_tag:
    fname = 'v0001_*'+str(tag)+'.dat'
    assert(len(glob.glob1(Base_folder,fname))==npart)

    fempty = True
    for i in range(npart):
        sname  = Base_folder + '/v0001_'+str(i)+'_'+str(tag)+'.dat'
        data   = np.loadtxt(sname, unpack=False)
        if(data != []):
            data = data.reshape(-1,7)
            if(fempty):
                fdata  = 1*data
                fempty = False
            else:
                fdata = np.concatenate((fdata,data))
    
    if(tag == surf_tag[0]):
        sfdata = fdata[fdata[:,0].argsort()]  
    else:
        sfdata = np.concatenate((sfdata,fdata[fdata[:,0].argsort()[::-1]]))       
    
cp  = -(sfdata[:,2:3] - p)/KE     
    
plt.figure()
plt.plot(sfdata[:,0],cp,'-o')
plt.xlim(0,1)
plt.xlabel('x')
plt.ylabel('-Cp')
plt.draw()
plt.savefig('Cp.pdf',dvi=300,bbox_inches='tight')
plt.clf()



