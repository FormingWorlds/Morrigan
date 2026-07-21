import numpy as np
from mass_loss import * 
import pdb 
from helper_functions import *

def merge_embryo(ap, Mp, Rp, Ms, ecc, v_c, live_status, b, atm_mass_fraction): #calculate orbital parameters post collision
    '''
    Function to compute mass, orbital separation, and eccentricity after a giant impact

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

    #atmospheric mass fraction lost from target body (here calculated from the merged body after collision)
    frac_lost = mass_loss(v_c, Mp[impactor], Mp_new, rho[impactor], rho[target], Rp[impactor], Rp[target], b)

    live_status[dead] = False #kill consumed planet
    ap[alive] = ap_new 
    Mp[alive] = Mp_new 
    ecc[alive] = ecc_new 
    atm_mass_fraction[alive] = Mp_new - (Mp_new * frac_lost) #update atmospheric mass fraction
    atm_mass_fraction[dead] = 0.0 

    return ap,Mp,ecc,live_status,atm_mass_fraction,frac_lost