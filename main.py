import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
from pylaplace import LaplaceCoefficient 

####INITIALISE SYSTEM######
N = 10 #number of planets
inner_edge = 0.005 #(AU)
outer_edge = 1 
M = 0.5 #planet mass (relative to Mearth)

def hill_sphere(a_i):
    return rH = a_i*(M/(3*(2*M)))**(1/3) #hill sphere for adjacent planets

a = np.array(N)
for i,r in enumerate(N): #allocate initial semi-major axes
    a_innermost = 0.1 #starting semi-major axis

    a[i] = a + 10*hill

    

a =  #semi-major axes (AU)
e = 0.01
I = 
#M,a,e,I allocated for each planet in 1,2..N
Ms = #stellar mass (relative to Msun)

max_time = 1e9 #evolution time (years)

#####CONSTANTS######
G = 6.67e-11 #m^3kg^-1s^-2


#####ANCILLARY FUNCTIONS#######

def del_ij(ei,ej,ai,aj)

def e_esc(Mi,Mj,Ri,Rj,aij)

def omega_min(eio,ai,ejo,aj)

def tau_col(Ri,Rj,v_esc,r_ran, eij)

def tau_scat(bij, aij, Ri, Rj, eij, e_esc)

def hi_dt(gj, betaj)

def ki_dt(gj, betaj)

####MAIN FUNCTIONS#########

def secular_ecc(h,k)
#big function for:
#1. calculating eccentricity via secular eccentricity procedure, calles pylaplace to calculate coefficients
#2. 

def tau_cross(M, nu, del_6, del_ov, P1)
#calculates time of next crossing (event) as well as the timescales of the events

def crossing_ecc(Mi,Mj, bij, ai, aj)

def event_path(tau_col, tau_scat)
#evaluates scattering vs collision, gives probability, p_col

def event_outcome(p_col)
#evaluates new mass, a, and e after either collision or scattering, 

def next_step()
#re-order planets based on eccentricities after an event, reset time to t = t+tau_cross + event time
#also evaluate if simulation is done (t = max_time or N = 2 stability criterion satisfied)