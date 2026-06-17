import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
from pylaplace import LaplaceCoefficient 

####INITIALISE SYSTEM######
N = 10 #number of planets
M = 1 #planet mass
a = 1 #semi-major axes
e = 
inner_edge = 0.005
Ms = #stellar mass

max_time = 1e9 #evolution time (years)

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

def tau_cross(M, nu, del_6, del_ov, P1)

def crossing_ecc(Mi,Mj, bij, ai, aj)

def event_path(tau_col, tau_scat)
#evaluates scattering vs collision, gives probability, p_col

def event_outcome(p_col)
#evaluates new mass, a, and e after either collision or scattering 

def next_step()
#re-order planets based on eccentricities after an event, reset time to t = t+tau_cross + event time