import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
from pylaplace import LaplaceCoefficient 
import pdb 

####INITIALISE SYSTEM######
N = 10 #number of planets
inner_edge = 0.005 #(AU)
outer_edge = 1 
M = 0.5 #planet mass (relative to Mearth)

def hill_sphere(a_i):
    return a_i*(M/(3*(2*M)))**(1/3) #hill sphere for adjacent planets

a = np.empty(N)
a[0] = 0.1
for i in range(N): #allocate initial semi-major axes
    pdb.set_trace
    a_previous = a[i-1] #starting semi-major axis
    hill = hill_sphere(a_previous)
    a[i] = a_previous + 10*hill

pdb.set_trace()

e = 0.01
#M,a,e,I allocated for each planet in 1,2..N
Ms = 1 #stellar mass (relative to Msun)

max_time = 1e9 #evolution time (years)

#####CONSTANTS######
G = 6.67e-11 #m^3kg^-1s^-2
