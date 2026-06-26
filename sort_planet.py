import numpy as np

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
    live_status = live_status[sort_order]

    return ap, Mp, ecc, Rp, live_status, interact, densities