import sys
import os
import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
from pylaplace import LaplaceCoefficient 
import pdb 

#####CONSTANTS######
G = 6.67e-11 #m^3kg^-1s^-2
M_sun = 1.9892e33 #grams
M_earth = 5.9736e27 #grams

####INITIALISE SYSTEM######
N = 10 #number of planets
inner_edge = 0.005 #(AU)
e = 0.01
Mp = 0.5 #planet mass (relative to Mearth)
Ms = 1*M_sun #stellar mass (relative to Msun)
max_time = 1e9 #evolution time (years)

def hill_sphere(a_i,M):
    return a_i * ((2*M) / (3 * Ms))**(1/3) #mutual hill radius for adjacent planets

def allocate_a(M):
    M *= 5.9736e27 #convert to Earth masses
    a = np.empty(N)
    a[0] = 0.1
    for i in range(1,N): #allocate initial semi-major axes
        a_previous = a[i-1] #starting semi-major axis
        hill = hill_sphere(a_previous,M)
        a[i] = a_previous + 10*hill #planets are spaced about by 10 hill radii
    return a

plt.scatter(allocate_a(0.1),range(N), label = 'M_p = 0.1')
plt.scatter(allocate_a(0.5), range(N), label = 'M_p = 0.5')
plt.scatter(allocate_a(1.0), range(N), label = 'M_p = 1.0')
plt.legend()
plt.xlabel('Initial orbital separation (AU)')
plt.ylabel('Planet number')
plt.tight_layout()
plt.savefig('initial_system.png',dpi = 300)
plt.close()

