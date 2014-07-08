


import numpy as np
import matplotlib.pyplot as plt
from numpy.lib.scimath import sqrt as csqrt


def TMM(period, section_length, ref_mod, delat, alpha,dz):
    
    T11 = np.cosh(alpha*dz)-1j*(delta/alpha)*np.sinh(alpha*dz)
    T22=np.cosh(alpha*dz)+1j*(delta/alpha)*np.sinh(alpha*dz)
    T12=-1j*(kapaac/alpha)*np.sinh(alpha*dz)
    T21=1j*(kapaac/alpha)*np.sinh(alpha*dz)

    return np.array([[T11,T12],[T21,T22]])



# constants - all dimensional untis in m

r_core = 5e-6 # radius of the core
lambda_cut = 700e-9 # wavelength in mm
n_core = 1.445 # refractive idnex of the core
n_clad = np.sqrt(-(2.405**2*lambda_cut**2)/(4*np.pi**2*r_core**2)+n_core**2) # claculate cladding oindex given cut off wavelenght and core radius and refractive index
d_n_max = 1e-4 # magnitude of index modulation


Bragg_wavelength = 750e-9
period = Bragg_wavelength/(2*n_core)
L = 2e-3
no_sections = 100
section_length = L/no_sections

lambda_min = 748.3e-9 
lambda_points = 500
lambda_range = 3e-9
lambda_max = lambda_min + lambda_range

wavelength = np.linspace(lambda_min,lambda_max,lambda_points)

ko = 2*np.pi/wavelength

#calculating effective index of fibre mode

vf=((2*np.pi/wavelength)*r_core*np.sqrt(n_core**2-n_clad**2))

u = (1+np.sqrt(2))*vf/(1+(4+vf**4)**0.25)
bfib = 1-u**2/vf**2

n_eff = np.sqrt(bfib*(n_core**2-n_clad**2)+n_clad**2)  
R=np.zeros((lambda_points))
dz = np.int(section_length/period)*period
Length = dz*no_sections

## flat top
#delta_n  = d_n_max*np.ones(no_sections)
#z = np.linspace(-Length/2,Length/2,no_sections)



## Guassian
#z = np.linspace(-Length/2,Length/2,no_sections)
#delta_n = d_n_max*np.exp(-(2*z/(0.4*Length))**2)

# Raised Cosine
z = np.linspace(-Length/2,Length/2,no_sections)
delta_n = 0.5*d_n_max*(1+np.cos(np.pi*(z/(0.1*Length))))



for i in range(0,lambda_points):
    MAT = np.eye(2)
    


    for m in range(no_sections-1,0,-1):
        kapadc = 4*np.pi*delta_n[m]/wavelength[i]
        kapaac = kapadc/2
        delta = kapadc+0.5*(2*ko[i]*n_eff[i]-2*np.pi/period)
        alpha = csqrt(kapaac**2-delta**2)

        T = TMM(period, dz, delta_n[m], delta, alpha,dz)
        
        MAT =np.dot(MAT,T)
        
    R[i] = np.abs(MAT[0,1]/MAT[0,0]) 

   
   
plt.figure(2)
plt.plot(wavelength,R)

plt.figure(3)
plt.plot(z,delta_n)
plt.show()


