from scipy.interpolate import griddata
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import numpy as np

#Shell environment
import os
import sys
import subprocess		# to send python variables to shell.

# Permutations
import itertools

def distance(_x):
    return np.sqrt(np.dot(_x,_x))

# Change here.
def wavefunction(_xp):
#     He atom.
    r1 = distance(_xp[0])
    r2 = distance(_xp[1])
    r12 = distance(_xp[0] - _xp[1])
    return np.exp( -2.0*(r1 + para[0]*r1*r1)/(1+r1) - 2.0*(r2+para[0]*r2*r2)/(1+r2) + 0.5*r12/(1+para[1]*r12) )
#     return np.exp( -para[0]*(r1 + r2) ) #Choice 2
    # SHO   - Harmonic oscillator
#     _sum = 0.0
#     for iP in range(len(_xp)):
#         _r = distance( _xp[iP] )
#         _sum = _sum + np.exp( -para[0]* _r**2 )
#     return _sum

#Change here.
def potential(_xp):
#     He atom                        
    r1 = distance(_xp[0])
    r2 = distance(_xp[1])
    r12 = distance(_xp[0] - _xp[1])
    return -2.0/r1 - 2.0/r2 + 1.0/r12
    #SHO
#     _sum = 0.0
#     for iP in range(len(_xp)):
#         _r = distance( _xp[iP] )
#         _sum = _sum + 0.5*_r**2
#     return _sum
    
def d2FdX(_xp):
    _dx = 0.0000001
    _sum = 0.0
    for iP in range(len(_xp)):
        for iD in range(len(_xp[0])):
            _xplus = np.copy(_xp)
            _xplus[iP][iD] = _xplus[iP][iD] + _dx
            _xminus = np.copy(_xp)
            _xminus[iP][iD] = _xminus[iP][iD] - _dx
            _sum = _sum +  ( wavefunction(_xplus) + wavefunction(_xminus) - 2*wavefunction(_xp) )/(_dx*_dx)
    return _sum

def localE(_xp):
    Vchi = potential(_xp) * wavefunction(_xp) # V(x) chi(x)
    Tchi = - 0.5*d2FdX(_xp) # T chi(x)
    return (Tchi + Vchi)*1.0/wavefunction(_xp) # (1/chi)*H*chi


#### Metropilis-Hasting symmetric sampling.
def vmcEnergy(_nVMC, _wfout):
    _Etotal = 0.0    
    _xsample = np.array( [ np.random.uniform(-1,1,nD) for i in range(nP)]) # Initial random starting point.
    for _ix in range(0, _nVMC):
        ########## MH sampling
        _xproposed = np.copy(_xsample)
        iP = np.random.randint(0, len(_xsample), 1) # Select a particle.
        _xproposed[iP] = _xsample[iP] + np.array( np.random.normal(0.0,1.0,len(_xsample[0]) ) )   # Generate new position around the old position using Gaussian.
        if np.random.uniform(0,1) <= (wavefunction(_xproposed)/wavefunction(_xsample))**2: # Accept/reject the new postition
            _xsample = np.copy(_xproposed)
        ###########
        # Local energy calculation
        #print(wavefunction(_xsample)**2, localE(_xsample))
        _Etotal = _Etotal + localE(_xsample)
        # Wavefunction out 
        if _wfout == 1:
            np.savetxt(f_wf, np.array( [ [distance(_xsample[0]), distance(_xsample[1]), distance(_xsample[0] -  _xsample[1]),  wavefunction(_xsample) ] ]), delimiter=" ", fmt="%0.12f")
    return _Etotal/_nVMC


#########
nP = 2 # Number of particles
nD = 3 # Number of dimension
p0 = [0.75, 3.3] # Initital parameter
p1 = [1.5, 8.0]
dpara = 0.1 

do_printwf = 1
do_scanpara = 0
#########
para = np.copy(p0)
########

######## Analyze wave function data.
f_wf = open('./wf.dat', 'w')   # To save wave fucntion data.
vmcEnergy(500000, do_printwf)
f_wf.close()
#######

############## Scan parameter range.
#result = np.array([]).reshape(0,len(para) + 1) # [para, E]

if do_scanpara == 1 :
    f_para = open('./vmcout.dat','w')
    while para[0] < p1[0]:
        para[1] = np.copy(p0[1])
        while para[1] < p1[1]:
            #e = vmcEnergy(500000, do_printwf)
            e = vmcEnergy(100000, do_printwf)
            print(para[0],para[1], e)
            #result = np.vstack( (result, np.array( [para[0], para[1], e ] ) ) )
            para[1] = para[1] + dpara 
            np.savetxt(f_para, np.array( [ [para[0], para[1], e ] ]) , delimiter=" ", fmt="%0.12f")
        para[0] = para[0] + dpara 
    f_para.close()
#################

# x0 = np.array( [ 0.0*np.random.uniform(-1,1,nD) for i in range(nP)])
# print("Here", x0, distance(x0[0]), d2FdX(x0), localE(x0))

### covergence test
# goodN = 200000
# dgoodN = 25000
# for iN in range(10):
#     print(goodN+iN*dgoodN, vmcEnergy(goodN+iN*dgoodN))
### 


