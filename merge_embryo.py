import numpy as np

def merge_embryo(ap, Mp, ecc, live_status): #calculate orbital parameters post collision
    Mp_new = sum(Mp) #eq 15
    ap_new = (Mp[0]*ap[0] + Mp[1]*ap[1])/Mp_new #eq 16

    cos_dvarpi = (ecc[0]**2*ap[0]**2 + ecc[1]**2*ap[1]**2 - (ap[0] - ap[1])**2)/(2*ecc[0]*ecc[1]*ap[0]*ap[1]) #eq 26
    min_dvarpi = np.arccos(cos_dvarpi) #eq 27
    dvarpi = np.random.uniform(min_dvarpi, 2*np.pi - min_dvarpi) #random range for dvarpi

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