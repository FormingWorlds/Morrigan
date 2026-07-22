"""
!!! info "`interaction_timescales.py`"
    Functions to calculate crossing timescale, viscous interaction, and collision timescale
    Author(s): Anna Grace Ulses
"""

import numpy as np 
from morrigan.constants import G
from morrigan.helper_functions import hill_sphere, kepler_period

def tau_cross_petit(a,Mp,Ms,ecc, N_affect): #evaluates every planetary triplet for instability
    '''
    Function to calculate timescale to instability as a result of 3-body mean motion resonances

    Parameters 
    ----------
    a : list 
        Semi-major axes of planetary triplet [m]
    Mp : list 
        Masses of planetary triplet [kg]
    Ms : float
        Stellar mass [kg]
    ecc : list
        Eccentricities of planetary triplet 
    N_affect : int
        Number of planets interacting 

    Returns
    -------
    tau_cross : float
        Instability timescale for planetary triplet [s]
    '''
    K = min(0.5*(N_affect-3) + 1, 3)
    alpha_01, alpha_12 = a[0]/a[1], a[1]/a[2]

    nu_01, nu_12 = alpha_01**1.5, alpha_12**1.5 #different from paper appendix B
    
    eta = (nu_01 * (1 - nu_12))/(1 - nu_01*nu_12)
    M = np.sqrt(Mp[0] * Mp[2] + Mp[1] * Mp[2] * eta**2 * alpha_01**(-2) + Mp[0]*Mp[1] * alpha_12**2 * (1 - eta)**2)/Ms
    delta_01, delta_12 = 1.0 - ecc[1] - (1.0 + ecc[0]) * alpha_01 , 1.0 - ecc[2] - (1.0 + ecc[1]) * alpha_12 #also different from paper, but noted in the code

    if delta_01 < 0.0 or delta_12 < 0.0:
        return 0.0

    delta = (delta_01 * delta_12) / (delta_01 + delta_12)
    delta_ov = (6.55 * K * M)**(1/4) * (eta*(1-eta))**(3/8)

    if delta >= delta_ov: #tau cross approaches infinity
        return 1e20

    log_arg = -np.log10((32 * np.sqrt(19) * M * np.sqrt(eta * (1 - eta)))/(3*np.sqrt(np.pi))) + np.log10(delta**6/(delta_ov**6 * (1 - (delta/delta_ov)**4))) + np.sqrt(-np.log(1 - (delta/delta_ov)**4))
    tau_cross = 10**(log_arg) * kepler_period(Mp[0],Ms,a[0])

    return tau_cross

def interaction_wrapper(ap, Mp, Ms, ecc, N_affect): #determine if system is stable, and if not, calculate timescale to instability
    '''
    Function to determine if triplet is stable, and if not, calculate timescale of instability

    Parameters
    ----------
    ap : list 
        Semi-major axes of triplet [m]
    Mp : list 
        Masses of triplet [kg]
    Ms : float
        Stellar mass [kg]
    ecc : list 
        Eccentricities of triplet 
    N_affect : int
        Number of planets interacting
    
    Returns: 
    1e20 if system is stable 
    Otherwise returns instability timescale from tau_cross_petit
    '''
    #N_affect is planets participating in crossing event
    aM = (Mp[0]*ap[0] + Mp[1]*ap[1]) / (Mp[0] + Mp[1])
    h = hill_sphere(aM,Mp[0]+Mp[1],Ms) / aM
    #stability criterion
    EJbef = 5.0/8.0*(ecc[0]**2 + ecc[1]**2)/h**2 - 3.0/8.0 * ((ap[0]-ap[1])/(h*aM))**2 + 4.5 #eq 28

    if N_affect <= 2 and EJbef < 0.0: #if system is stable return 'infinite' stability, 2 planets can still interact
        return 1e20
        
    #otherwise return the crossing timescale from Petit 2020
    return tau_cross_petit(ap, Mp, Ms, ecc, N_affect)

def tau_vis(ap,Mp,Rp,Ms,ecc): #viscous relaxation timescale for an interacting planetary pair
    '''
    Viscous relaxation timescale for an interacting planetary pair

    Parameters
    ----------
    ap : list
        Semi-major axes of interacting pair [m]
    Mp : list 
        Masses of interacting pair [kg]
    Rp : list 
        Radii of interacting pair [m]
    Ms : float
        Stellar mass [kg]
    ecc : list
        Eccentricities of interacting pair 

    Returns
    -------
    timescale : float 
        Viscous relaxation timescale (when system settles down after scattering)
    '''
    mu_a = sum(ap)/2 #average semi-major axis of interacting pair
    mu_e = np.sqrt(ecc[0]**2 + ecc[1]**2)
    M_T = sum(Mp) #sum of masses
    impact_parameter = abs(ap[1] - ap[0]) 

    #eccentricities at the onset of crossing
    ecross_i = (np.sqrt(Mp[1]) * impact_parameter)/((np.sqrt(Mp[1]) * ap[0]) + np.sqrt(Mp[0]) * ap[1]) #eq 6
    ecross_j = (np.sqrt(Mp[0]) * impact_parameter)/((np.sqrt(Mp[0]) * ap[0]) + np.sqrt(Mp[1]) * ap[1]) #implied eq 6
    ecross = [ecross_i, ecross_j]

    rep_e = 0.5 * max((sum(ecross)), sum(ecc)) #eq 23 used to calculate lambda in eq 12
    kep_vel = np.sqrt((G * Ms)/mu_a) #eq 10
    ran_vel = rep_e * kep_vel #before eq 9
    n = 1/(2 * np.pi * rep_e * mu_a**2 * impact_parameter) #eq 8

    timescale = ran_vel**3 / (n * np.pi * G**2 * M_T**2) # matches Fortran (no factor of 3)
    return timescale

def tau_col(ap,Mp,Rp,Ms,ecc):
    '''
    Function to calculate duration of a collision event

    Parameters
    ----------
    ap : list
        Semi-major axes of interacting pair [m]
    Mp : list 
        Masses of interacting pair [kg]
    Rp : list 
        Radii of interacting pair [m]
    Ms : float
        Stellar mass [kg]
    ecc : list
        Eccentricities of interacting pair 

    Returns
    -------
    timescale : float
        Duration of collisional event [s]
    '''
    mu_a = sum(ap)/2 #average semi-major axis of interacting pair
    mu_e = np.sqrt(ecc[0]**2 + ecc[1]**2)
    M_T = sum(Mp) #sum of masses
    R_T = sum(Rp) #sum of radii
    impact_parameter = abs(ap[1] - ap[0])

    #eccentricities at the onset of crossing
    ecross_i = (np.sqrt(Mp[1]) * impact_parameter)/((np.sqrt(Mp[1]) * ap[0]) + np.sqrt(Mp[0]) * ap[1]) #eq 6
    ecross_j = (np.sqrt(Mp[0]) * impact_parameter)/((np.sqrt(Mp[0]) * ap[0]) + np.sqrt(Mp[1]) * ap[1]) #implied eq 6
    ecross = [ecross_i, ecross_j]

    rep_e = 0.5 * max(sum(ecross), sum(ecc))
    kep_vel = np.sqrt((G * Ms) / mu_a)
    ran_vel = rep_e * kep_vel
    esc_vel = np.sqrt(2.0 * G * M_T / R_T)
    n = 1.0 / (2.0 * np.pi * rep_e * mu_a**2 * impact_parameter)

    #only difference from tau_vis here is the timescale
    timescale = n * np.pi * (R_T)**2 * (1 + esc_vel**2/ran_vel**2) * ran_vel #eq 11

    return 1/timescale