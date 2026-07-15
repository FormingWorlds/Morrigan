import numpy as np
from mass_loss import * 
import pdb 

def merge_embryo(ap, Mp, Ms, ecc, v_c, live_status): #calculate orbital parameters post collision
    #Mp is a list containing interacting pair masses
    Mp_new = sum(Mp) #eq 15
    ap_new = Mp_new/(Mp[0]/ap[0] + Mp[1]/ap[1]) #eq 16

    if ap[1]*(1.0 - ecc[1]**2) < ap[0]*(1.0 - ecc[0]**2):
        min_dvarpi = 0.0
    else:
        cosdvarpi = ((ecc[0]*ap[0])**2 + (ecc[1]*ap[1])**2 - (ap[1]-ap[0])**2) / (2.0 * ecc[0] * ecc[1] * ap[0] * ap[1])
        #cosdvarpi = min(1.0, cosdvarpi)
        cosdvarpi = max(-1.0, min(1.0, cosdvarpi))
        min_dvarpi = np.arccos(cosdvarpi)

    dvarpi = np.random.uniform(min_dvarpi, 2.0 * np.pi - min_dvarpi)
    ecc_new = np.sqrt(((Mp[0]**2*ecc[0]**2) + (Mp[1]**2*ecc[1]**2) + 2*Mp[0]*Mp[1]*ecc[0]*ecc[1]*np.cos(dvarpi)) / Mp_new**2) #eq 17

    if Mp[0] >= Mp[1]: #larger planet consumes smaller one 
        alive, dead = 0,1
    else:
        alive, dead = 1,0

    live_status[dead] = False 
    ap[alive] = ap_new 
    Mp[alive] = Mp_new 
    ecc[alive] = ecc_new 

    return ap,Mp,ecc,live_status