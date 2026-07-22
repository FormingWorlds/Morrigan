"""
!!! info "`merge_embryo.py`"
    Handles merging events. 
    collision_velocity: calculates collision velocity between target and impactor
    merge_embryo: Computes resulting mass, orbital separation, eccentricity, and atmospheric mass loss
    Author(s): Anna Grace Ulses
"""

import numpy as np
from mass_loss import * 
import pdb 
from helper_functions import *

def collision_velocity(ap, Mp, Rp, Ms, ecc):
    '''
    Velocity of collision between target and impactor during merge events

    Parameters:
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
    v_c : float
        Collision velocity between target and impactor during merging [km/s]
    '''
    mu_a = sum(ap) / 2
    kep_vel = np.sqrt(G * Ms / mu_a)
    rep_e = np.sqrt(ecc[0]**2 + ecc[1]**2)  # same as eij in orbit_cross_K25
    v_inf = rep_e * kep_vel
    v_esc = np.sqrt(2 * G * (Mp[0] + Mp[1]) / (Rp[0] + Rp[1]))
    v_c = np.sqrt(v_inf**2 + v_esc**2)
    return v_c


def merge_embryo(ap, Mp, Rp, Ms, ecc, v_c, live_status, b, atm_mass_fraction): #calculate orbital parameters post collision
    '''
    Function to compute mass, orbital separation, eccentricity, and fractional atmospheric loss after a giant impact

    Parameters 
    ----------
    ap : list
        Semi major axes of colliding pair [m]
    Mp : list 
        Masses of colliding pair [kg]
    Rp : list 
        Radii of colliding pair [m]
    ecc : list 
        Secular eccentricity of colliding pair just before collision
    v_c : float 
        Collision velocity [m/s]
    live_status : list 
        Live_status of colliding pair 
    b : float 
        Impact parameter, defined as sin(beta) where beta is impact angle in [deg]
    atm_mass_fraction : list 
        Atmospheric mass fractions of interacting pair
    
    Returns
    -------
    ap : float 
        Updated semi-major axis of surviving planet [m]
    Mp : float
        Updated mass of surviving planet [kg]
    ecc : float 
        Updated eccentricity of surviving planet 
    live_status : list
        Sets consumed planet's status to 'False'
    atm_mass_fraction: float 
        Remaining atmospheric mass fraction of target after collision
    frac_lost : float
        Fraction of surviving planets atmosphere that was lost during collision

    '''
    #Mp is a list containing interacting pair masses
    Mp_new = sum(Mp) #eq 15
    ap_new = Mp_new/(Mp[0]/ap[0] + Mp[1]/ap[1]) #eq 16

    if ap[1]*(1.0 - ecc[1]**2) < ap[0]*(1.0 - ecc[0]**2):
        min_dvarpi = 0.0
    else:
        cosdvarpi = ((ecc[0]*ap[0])**2 + (ecc[1]*ap[1])**2 - (ap[1]-ap[0])**2) / (2.0 * ecc[0] * ecc[1] * ap[0] * ap[1])
        cosdvarpi = max(-1.0, min(1.0, cosdvarpi))
        min_dvarpi = np.arccos(cosdvarpi)

    dvarpi = np.random.uniform(min_dvarpi, 2.0 * np.pi - min_dvarpi)
    ecc_new = np.sqrt(((Mp[0]**2*ecc[0]**2) + (Mp[1]**2*ecc[1]**2) + 2*Mp[0]*Mp[1]*ecc[0]*ecc[1]*np.cos(dvarpi)) / Mp_new**2) #eq 17

    if Mp[0] >= Mp[1]: #larger planet consumes smaller one 
        alive, dead = 0,1
    else:
        alive, dead = 1,0

    target, impactor = alive, dead #by convention the surviving (larger) body is the target

    rho = Mp / ((4.0/3.0) * np.pi * Rp**3) #density

    #pre-merge atmosphere masses [kg], derived from the CARRIED-OVER fraction of each body
    #(this is what actually threads atmosphere state through successive collisions -- a planet
    #that already lost atmosphere in an earlier merger has a smaller atm_mass_fraction here)
    atm_mass_target_before = Mp[target] * atm_mass_fraction[target]
    atm_mass_impactor_before = Mp[impactor] * atm_mass_fraction[impactor]
 
    #impactor's atmosphere is added to the target's before mass_loss() is applied
    atm_mass_combined_before = atm_mass_target_before + atm_mass_impactor_before

    #fraction of that COMBINED atmosphere lost in this collision (Kegerreis et al. 2020 scaling)
    frac_lost = mass_loss(v_c, Mp[impactor], Mp[target], rho[impactor], rho[target], Rp[impactor], Rp[target], b)
    frac_lost = min(max(frac_lost, 0.0), 1.0) #the fitted scaling is only defined/meaningful on [0,1]
 
    atm_mass_after = atm_mass_combined_before * (1.0 - frac_lost)
    atm_mass_lost = atm_mass_combined_before - atm_mass_after
 
    Mp_new = Mp_new - atm_mass_lost #reduce the merged mass by the atmosphere actually lost (rock is conserved)
 
    live_status[dead] = False #kill consumed planet
    ap[alive] = ap_new 
    Mp[alive] = Mp_new 
    ecc[alive] = ecc_new 
    atm_mass_fraction[alive] = atm_mass_after / Mp_new #store back as a FRACTION of the new total mass
    atm_mass_fraction[dead] = 0.0 

    return ap,Mp,ecc,live_status,atm_mass_fraction,frac_lost