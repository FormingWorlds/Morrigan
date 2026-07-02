import numpy as np 
from constants import *

def kepler_period(Mp, Ms, a):
    P_squared = (4*np.pi**2*a**3)/(G*Ms)   # matches Fortran: only stellar mass
    return np.sqrt(P_squared)

def rayleigh(sigma, xmin):
    Umin = 1.0 - np.exp(-0.5 * xmin**2 / sigma**2)
    dum = np.random.uniform(Umin, 1.0 - 1e-10)
    return sigma * np.sqrt(-2.0 * np.log(1.0 - dum))

def esc_ecc(Ms,M1,M2,R1,R2,a):
    num = np.sqrt(2*G*(M1+M2)/(R1+R2))
    denom = np.sqrt((G*Ms)/a)
    return num/denom

def planet_radius(mass,density):
    return ((3*mass)/(4 * np.pi*density))**(1/3)

def hill_sphere(a,M,Ms):
    return a * ((M) / (3 * Ms))**(1/3) #mutual hill radius for adjacent planets