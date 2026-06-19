import sys
import os
import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
from pylaplace import LaplaceCoefficient 
import pdb 
from astropy.table import Table
from astropy.io import ascii

#####CONSTANTS######
G = 6.67e-11 #m^3kg^-1s^-2
M_sun = 1.9892e30 #kg
M_earth = 5.9736e24 #kg

####INITIAL PARAMETERS######
N = 10 #number of planets
e = 0.01
Mp = 0.5*M_earth #planet mass (relative to Mearth)
Ms = 1*M_sun #stellar mass (relative to Msun)
max_time = 1e9*365*24*60*60 #evolution time (seconds)
rho_p = 5500 #planet density kg/m^3

####ALLOCATE PARAMETERS FOR THE SYSTEM###
def hill_sphere(a_i,M):
    return a_i * ((2*M) / (3 * Ms))**(1/3) #mutual hill radius for adjacent planets

def allocate_a(M):
    a = np.empty(N)
    a[0] = 0.1 #AU
    for i in range(1,N): #allocate initial semi-major axes
        a_previous = a[i-1] #starting semi-major axis
        hill = hill_sphere(a_previous,M)
        a[i] = a_previous + 10*hill #planets are spaced out by 10 hill radii
    return a*1.5e11 #convert to [m] to stay in SI!

def planet_radius(mass,density):
    return ((3*mass)/(4 * np.pi*density))**(1/3)

a = allocate_a(Mp)
ecc = np.full(N,e)
densities = np.full(N, rho_p)
masses = np.full(N, Mp)
live_status = np.ones(N, dtype = bool) #set initial status of planets, all are live by definition at the start
Rp = [planet_radius(i, j) for i,j in zip(masses,densities)]

parameter_names = ['a_AU','e','Mp','Rp','live_status']
system_information = Table([a/1.5e11,ecc,masses,Rp,live_status],names = parameter_names)
ascii.write(system_information, 'initial_system.csv', format = 'fixed_width', overwrite = True)

######SECULAR ECCENTRICITY SECTION#####

n = np.empty(N)
for i in range(N):
    n[i] = np.sqrt(G*Ms/a[i]**3)

A = np.zeros((N,N)) #empty interaction matrix
laplace = LaplaceCoefficient(method = 'Brute')
for i in range(N): 
    for j in range(N): 
        if i == j:
            continue #skip self-interactions
        if a[i] < a[j]:
            alpha = a[i]/a[j]
            alpha_bar = alpha 
        else:
            alpha = a[j]/a[i]
            alpha_bar = 1

        #calculate laplacian coefficients, b1,3/2 and b2,3/2
        #(a,s,m,p,q)
        coeff_m1 = laplace(alpha, 3/2, 1, 1, 1) #m = 1 for A_ii
        coeff_m2 = laplace(alpha, 3/2, 2, 1, 1) #m = 2 for A_ij

        factor = n[i] * 0.25 * masses[j]/(Ms + masses[i]) * alpha * alpha_bar
        A[i,i] += factor * coeff_m1
        A[i,j] = -factor * coeff_m2

# Set initial angles randomly
varpi = np.random.uniform(0.0, 2.0 * np.pi, N)

h0 = ecc * np.sin(varpi)
k0 = ecc * np.cos(varpi)

# Solve the eigenvalue problem
g, S = np.linalg.eig(A)
g = np.real(g)
S = np.real(S)

# Solve for integration constants using S
Csinb = np.linalg.solve(S, h0)
Ccosb = np.linalg.solve(S, k0)

# Calculate amplitudes (C) and phase angles (beta)
C = np.sqrt(Csinb**2 + Ccosb**2)
beta = np.arctan2(Csinb, Ccosb)

# Scale eigenvectors by the amplitudes (columns of ecc_vec)
ecc_vec = S * C

#####ANCILLARY FUNCTIONS##########

def kepler_P(Mp,a) #period of planetary orbit, used to calculate tau_cross
    P_squared = (4*np.pi**2*a**3)/(G*(Mp+Ms))
    return np.sqrt(P_squared)

def esc_ecc(M1,M2,R1,R2,a):
    numerator = np.sqrt(2*G*(M1+M2)/(R1+R2))
    denom = np.sqrt((G*Ms)/s)
    return numerator/denominator

def tau_cross_petit: #evaluates for a planetary triplet?
    K = 
    alpha_01, alpha_12 = a[0]/a[1], a[1]/a[2]
    nu_01, nu_12 = kepler_period(M[0],a[0])/kepler_period(M[1],a[1]) , kepler_period(M[1],a[1])/kepler_period(M[2],a[2])
    eta = (nu_01 * (1 = nu_12))/(1 - nu_01*nu_12)
    M = np.sqrt(Mp[0] * Mp[2] + Mp[1] * Mp[2] * eta**2 * alpha_01**(-2) + Mp[0]*Mp[2] * alpha_12**2 * (1 - eta)**2)/Ms
    delta_01, delta_12 = ((1-ecc[1])*a[1] - (1+ecc[0])*a[0])/a[1] , ((1-ecc[1])*a[1] - (1+ecc[2])*a[2])/a[1]
    delta = (delta_01 * delta_12) / (delta_01 + delta_12)
    delta_ov = (6.55 * K * M)**(1/4) (eta*(1-eta))**(3/8)
    delta_6 = 

    log_arg = -np.log10((32 * np.sqrt(19) * M * np.sqrt(eta * (1 - eta)))/(3*np.sqrt(np.pi)))\
    + np.log10(delta**6/(delta_ov**6 * (1 - (delta/delta_ov)**4))) + np.sqrt(-np.log(1 - (delta/delta_ov)**4))

    tau_cross = np.exp(log_arg) * kepler_period(M[0],a[0])

    return tau_cross

