"""
!!! info "`sort_planet.py`"
    Function to clean up planetary systems and sort by semi-major axis after every event
    Author(s): Anna Grace Ulses
"""

import numpy as np

def sort_planet(ap, Mp, ecc, Rp, live_status, interact, densities, planet_id): #clean up system after event
    '''
    Removes dead planets and sorts system by increasing semi-major axis 

    Parameters
    ----------
    ap : list
        Semi-major axes of all planets [m]
    Mp : list 
        Masses of all planets [m]
    ecc : list
        Eccentricities of all planets 
    Rp : list 
        Radii of all planets [m]
    live_status: bool list 
        Whether a planet is surviving (True), or was ejected or consumed (False)
    interact: bool list
        Indices of interacting planets at this timestep 
    densities: list
        Densities of all planets [kg/m^3]
    planet_id:
        Persistent planet ID to track evolutions
    
    Returns
    -------
    Sorted arrays for ap, Mp, ecc, Rp, atm_mass_fraction, live_status, interact, densities, and planet_id

    '''
    #remove ejected or dead ones 
    live_planets = live_status
    ap = ap[live_planets]
    Mp = Mp[live_planets]
    Rp = Rp[live_planets]
    ecc = ecc[live_planets]
    densities = densities[live_planets]
    planet_id = planet_id[live_planets]

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
    planet_id = planet_id[sort_order]

    return ap, Mp, ecc, Rp, live_status, interact, densities, planet_id