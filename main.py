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
np.random.seed(1) #for reproducability

####ALLOCATE PARAMETERS FOR THE SYSTEM###
def hill_sphere(a_i,M):
    return a_i * ((M) / (3 * Ms))**(1/3) #mutual hill radius for adjacent planets

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

#actually initialising system here with arrays for every parameter
a = allocate_a(Mp)
ecc = np.full(N,e)
densities = np.full(N, rho_p)
masses = np.full(N, Mp)
live_status = np.ones(N, dtype = bool) #set initial status of planets, all are live by definition at the start
interact = np.ones(N, dtype = bool) #stores the indices of which planets are participating in an event
Rp = [planet_radius(i, j) for i,j in zip(masses,densities)]

parameter_names = ['a_AU','e','Mp','Rp','live_status']
system_information = Table([a/1.5e11,ecc,masses,Rp,live_status],names = parameter_names)
ascii.write(system_information, 'initial_system.csv', format = 'fixed_width', overwrite = True) #store initial system information

######SECULAR ECCENTRICITY SECTION#####


def secular_solution(ap, Mp, ecc, Rp, N):
    varpi = np.random.uniform(0.0, 2.0 * np.pi, N)

    h0 = ecc * np.sin(varpi)
    k0 = ecc * np.cos(varpi)

    mean_motion = np.sqrt(G*Ms/ap**3)
    A = np.zeros((N,N)) #empty interaction matrix

    laplace = LaplaceCoefficient(method = 'Brute') #to calculate laplace coefficients
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

    #solve the eigenvalue problem
    g, S = np.linalg.eig(A)
    g = np.real(g)
    S = np.real(S)

    #solve for integration constants using S
    Csinb = np.linalg.solve(S, h0)
    Ccosb = np.linalg.solve(S, k0)

    #calculate amplitudes (C) and phase angles (beta)
    C = np.sqrt(Csinb**2 + Ccosb**2)
    beta = np.arctan2(Csinb, Ccosb)

    #scale eigenvectors by the amplitudes (columns of ecc_vec)
    ecc_vec = S * C

    return ecc_vec, g, beta

#####HELPER FUNCTIONS##########

def hill_sphere_mutual(M_sum, a_mean):
    return a_mean * (M_sum / (3.0 * Ms))**(1/3)

def kepler_P(Mp,a) #period of planetary orbit, used to calculate tau_cross
    P_squared = (4*np.pi**2*a**3)/(G*(Mp+Ms))
    return np.sqrt(P_squared)

def esc_ecc(M1,M2,R1,R2,a):
    num = np.sqrt(2*G*(M1+M2)/(R1+R2))
    denom = np.sqrt((G*Ms)/s)
    return num/denom

def tau_cross_petit(a,Mp,ecc, N_affect): #evaluates every planetary triplet for instability
    K = min(0.5*(N-3) + 1, 3)
    alpha_01, alpha_12 = a[0]/a[1], a[1]/a[2]
    nu_01, nu_12 = kepler_period(M[0],a[0])/kepler_period(M[1],a[1]) , kepler_period(M[1],a[1])/kepler_period(M[2],a[2])
    eta = (nu_01 * (1 = nu_12))/(1 - nu_01*nu_12)
    M = np.sqrt(Mp[0] * Mp[2] + Mp[1] * Mp[2] * eta**2 * alpha_01**(-2) + Mp[0]*Mp[2] * alpha_12**2 * (1 - eta)**2)/Ms
    delta_01, delta_12 = ((1-ecc[1])*a[1] - (1+ecc[0])*a[0])/a[1] , ((1-ecc[1])*a[1] - (1+ecc[2])*a[2])/a[1]
    delta = (delta_01 * delta_12) / (delta_01 + delta_12)
    delta_ov = (6.55 * K * M)**(1/4) (eta*(1-eta))**(3/8)

    log_arg = -np.log10((32 * np.sqrt(19) * M * np.sqrt(eta * (1 - eta)))/(3*np.sqrt(np.pi))) + np.log10(delta**6/(delta_ov**6 * (1 - (delta/delta_ov)**4))) + np.sqrt(-np.log(1 - (delta/delta_ov)**4))
    tau_cross = np.exp(log_arg) * kepler_period(M[0],a[0])

    return tau_cross

def interaction_wrapper(ap, Mp, ecc, N_affect): #determine if system is stable, and if not, calculate timescale to instability
    #N_affect is planets participating in crossing event
    aM = (Mp[0]*ap[0] + Mp[1]*ap[1]) / (Mp[0] + Mp[1])
    h = hill_sphere_mutual(Mp[0]+Mp[1], aM) / aM
    #stability criterion
    EJbef = 5.0/8.0*(ecc[0]**2 + ecc[1]**2)/h**2 - 3.0/8.0 * ((ap[0]-ap[1])/(h*aM))**2 + 4.5 #eq 28

    if N_affect <= 2 and EJbef < 0.0: #if system is stable return 'infinite' stability
        return 1e20
        
    #otherwise return the crossing timescale from Petit 2020
    return tau_cross_Petit(ap, Mp, ecc, N_affect)

def tau_vis(ap,Mp,Rp,ecc): #viscous relaxation timescale for an interacting planetary pair
    mu_a = sum(ap)/2 #average semi-major axis of interacting pair
    mu_e = sum(ecc)/2
    M_T = sum(Mp) #sum of masses
    impact_parameter = abs(ap[1] - ap[0]) 

    #eccentricities at the onset of crossing
    ecross_i = (np.sqrt(M[1]) * impact_parameter)/((np.sqrt(M[1]) * a[0]) + np.sqrt(M[0]) * a[1]) #eq 6
    ecross_j = (np.sqrt(M[0]) * impact_parameter)/((np.sqrt(M[0]) * a[0]) + np.sqrt(M[1]) * a[0]) #implied eq 6
    ecross = [ecross_i, ecross_j]

    rep_e = max((sum(ecross)), sum(ecc)) #eq 23 used to calculate lambda in eq 12
    kep_vel = np.sqrt((G * Ms)/mu_a)
    ran_vel = rep_e * kep_vel
    n = 1/(2 * np.pi * rep_e * mu_a**2 * impact_parameter)

    timescale = n * np.pi * G**2 * ran_vel * 3
    return 1/timescale

def tau_col(ap,Mp,Rp,ecc):
    mu_a = sum(ap)/2 #average semi-major axis of interacting pair
    mu_e = sum(ecc)/2
    M_T = sum(Mp) #sum of masses
    R_T = sum(Rp) #sum of radii
    impact_parameter = abs(ap[1] - ap[0])

    #eccentricities at the onset of crossing
    ecross_i = (np.sqrt(M[1]) * impact_parameter)/((np.sqrt(M[1]) * a[0]) + np.sqrt(M[0]) * a[1]) #eq 6
    ecross_j = (np.sqrt(M[0]) * impact_parameter)/((np.sqrt(M[0]) * a[0]) + np.sqrt(M[1]) * a[0]) #implied eq 6
    ecross = [ecross_i, ecross_j]

    rep_e = 0.5 * max(sum(ecc_cross), sum(ecc))
    kep_vel = np.sqrt(G * Ms / map_val)
    ran_vel = mecc * vK
    esc_vel = np.sqrt(2.0 * G * sMp / sRp)
    n = 1.0 / (2.0 * np.pi * mecc * map_val**2 * bimp)

    #only difference from tau_vis here is the timescale
    timescale = n * np.pi * (R_T)**2 * (1 + esc_vel**2/ran_vel**2) * ran_vel

    return 1/timescale

def crossing_pair(ap, Mp, Rp, ecc, ecc_vec, g, beta, interact, N, t, t_ref): #identify crossing pair from triplet, return pair and t_event

def merge_embryo(ap, Mp, ecc, live_status): #calculate orbital parameters post collision
    Mp_new = sum(Mp) #eq 15
    ap_new = (Mp[0]*ap[0] + Mp[1]*ap[1])/Mp_new #eq 16

    cos_dvarpi = (ecc[0]**2*ap[0]**2 + ecc[1]**2*ap[1]**2 - (ap[0] - ap[1])**2)/(2*ecc[0]*ecc[1]*ap[0]*ap[1]) #eq 26
    min_dvarpi = np.arccos(cos_dvarpi) #eq 27
    dvarpi = np.random.uniform(min_dvarpi, 2*np.pi - min_dvarpi) #random range for dvarpi

    ecc_new = np.sqrt(((Mp[0]**2*ecc[0]**2) + (Mp[1]**2*ecc[1]**2) + 2*Mp[0]*Mp[1]*ecc[0]*ecc[1]*np.cos(dvarpi)) / M_new**2) #eq 17

    if Mp[0] >= Mp[1]: #larger planet consumes smaller one 
        alive, dead = 0,1
    else:
        alive, dead = 1,0

    live_status[dead] = False 
    ap[alive] = ap_new 
    Mp[alive] = Mp_new 
    ecc[alive] = ecc_new 

    return ap,Mp,ecc,live_status

def orbit_cross_K25(ap, Mp, Rp, ecc, interact, live_status, N, icross): #determine outcome of crossing event
    #now working with an interacting pair of planets i,j
    jcross = icross + 1 #sets indices of interacting pair 
    mean_ap = 0.5*(ap[icross] + ap[jcross]) #average semi-major 
    e_esc = esc_ecc(Mp[icross],Mp[jcross],Rp[icross],Rp[jcross],mean_ap) #escape eccentricity

    ecross_i = (ap[jcross] - ap[icross]) * np.sqrt(Mp[jcross]) / (ap[icross] * np.sqrt(Mp[jcross]) + ap[jcross] * np.sqrt(Mp[icross])) #eq 6 again
    ecross_j = (ap[jcross] - ap[icross]) * np.sqrt(Mp[icross]) / (ap[icross] * np.sqrt(Mp[jcross]) + ap[jcross] * np.sqrt(Mp[icross]))

    ecc_cross = [ecross_i, ecross_j] #store together 

    ecc_dum = [max(ecc_cross[0], ecc[icross]), max(ecc_cross[1], ecc[jcross])] #eq 23


    #calculates collosion probability, Pcol
    #calles merge_embryo if a collision happens
    #in case of scattering excites the ecc and a_ij, if an ecc > 1, planet is ejected 

def sort_planet(ap, Mp, ecc, Rp, live_status, interact, densities): #clean up system after event
    #remove ejected or dead ones 
    live_planets = live_status
    ap = ap[live_planets]
    Mp = Mp[live_planets]
    Rp = Rp[live_planets]
    ecc = ecc[live_planets]
    densities = densities[live_planets]

    interact = interact[live_planets]
    live_status = live_status[live_planets]

    #sort by new semi-major axes
    sort_order = np.argsort(ap)
    ap = ap[sort_order]
    Mp = Mp[sort_order]
    ecc = ecc[sort_order]
    Rp = Rp[sort_order]
    interact = interact[sort_order]
    densities = densities[sort_order]
    live_status = live[sort_order]

    return ap, Mp, ecc, Rp, live_status, interact, densities
