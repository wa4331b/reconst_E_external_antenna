import numpy as np
import scipy.constants

ur  =1.     #permeability
er  =70.6   #permittivity of dielectric
f   =402e6  #frequency
lamda    =scipy.constants.c/f    #wavelength in vacuum
lamdag   =lamda/np.sqrt(er)      #wavelength in dielectric
h_lamdag =0.2   #height of coil normalized by lamdag
d_lamdag =0.108 #diameter of coil normalized by lamdag
t        =0.55  #tapering coefficient
r        =0.0006 #radius of wire
N        =7    #no of turns
#L        =np.pi*np.power(r,2)*N  #length of wire
L        =np.pi*np.power(r,2)  #length of wire
s        =58e6                   #conductivity wire

rd  =np.sqrt(ur/er)*20*np.power(np.pi,2)*np.power((h_lamdag),2) #small dipole resistance
rl  =20*np.power(np.pi,6)*np.sqrt(ur/er)*np.power((d_lamdag),4) #small loop resistance
ro  =t*(L/r)*np.sqrt(30/s*lamdag)*np.power(ur/er,0.25)     #ohmic resistance

print ('small dipole resistance is', rd)
print ('small loop resistance is', rl)
print ('ohmic resistance is', ro)