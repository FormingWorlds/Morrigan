"""
!!! info "`orbit_cross_K25.py`"
    Function that modifies arrays of orbital separation, mass, radius, eccentricity, interaction status, and live status of planets based on whether a collision or scattering (or ejection) event happens
    Author(s): Anna Grace Ulses
"""

import numpy as np

from morrigan.helper_functions import esc_ecc, hill_sphere, rayleigh
from morrigan.merge_embryo import collision_velocity, merge_embryo
 
def orbit_cross_K25(ap, Mp, Rp, Ms, atm_mass_fraction, impact_parameter, ecc, interact, live_status, N, planet_id, icross): #determine outcome of crossing event
    '''
    Function to determine what happens once an orbit crossing event has been established

    Parameters
    ----------
    ap : list 
        Semi major axes of interacting pair [m]
    Mp : list 
        Masses of interacting pair [kg]
    Rp : list 
        Radii of interacting pair [m]
    Ms : float
        Stellar mass [kg]
    atm_mass_fraction : list 
        Atmospheric mass fractions of interacting pair 
    impact_parameter: float
        Defined as sin(impact angle). Describes angle of contact between target and impactor
    ecc : list
        Eccentricities of interacing pair 
    interact : bool
        Interaction status of interacting pair (set to True)
    live_status : bool
        Live status of interacting pair (starts as True and will be modified depending on scattering vs. merge events)
    N : int
        Number of planets 
    planet_id : list 
        Persistent ID number for interacting planets to track evolutions
    icross: int
        Index of inner planet in interacting pair
    
    Returns
    -------
    Modified arrays for ap, Mp, Rp, atm_mass_fraction, ecc, interact, live_status, N following scattering or merging 
    merge_record : dict
        Stores information about targets and impactors for each merging event
    '''

    #now working with an interacting pair of planets i,j
    #modifies the arrays of ap,Mp,Rp,ecc,interact,live_status based on what happens
    jcross = icross + 1 #sets indices of interacting pair, aj>ai always, +1 to be able to index a pair later
 
    #mass-weighted mean semi-major axis of the pair (used for EJbef stability condition and the ejection formula)
    aM = (Mp[icross]*ap[icross] + Mp[jcross]*ap[jcross]) / (Mp[icross] + Mp[jcross])
    h = hill_sphere(aM, Mp[icross]+Mp[jcross], Ms) / aM
    EJbef = 5.0/8.0*(ecc[icross]**2 + ecc[jcross]**2)/h**2 - 3.0/8.0 * ((ap[icross]-ap[jcross])/(h*aM))**2 + 4.5
 
    #2-body case: if the pair is Jacobi-stable and not currently overlapping, nothing happens
    if N == 2 and (EJbef < 0.0 and (1.0+ecc[icross])*ap[icross] < (1.0-ecc[jcross])*ap[jcross]):
        return
 
    mean_ap = (ap[icross] + ap[jcross])/2 #simple average semi-major axis, used for e_esc
    e_esc = esc_ecc(Ms,Mp[icross],Mp[jcross],Rp[icross],Rp[jcross],mean_ap) #escape eccentricity
 
    ecross_i = (ap[jcross] - ap[icross]) * np.sqrt(Mp[jcross]) / (ap[icross] * np.sqrt(Mp[jcross]) + ap[jcross] * np.sqrt(Mp[icross])) #eq 6 again
    ecross_j = (ap[jcross] - ap[icross]) * np.sqrt(Mp[icross]) / (ap[icross] * np.sqrt(Mp[jcross]) + ap[jcross] * np.sqrt(Mp[icross]))
 
    ecc_cross = [ecross_i, ecross_j] #store together 
    #minimum eccentricity of interacting planets to cause a collision
 
    ecc_encounter = [max(ecc_cross[0], ecc[icross]), max(ecc_cross[1], ecc[jcross])] #eq 23
    #enforces that the planets ACTUALLY interact (actual eccentricites from ecc[] are secular, not at exact crossing time)
    #crossing_pair is still just a statistical prediction
    #ensures the eccentricities used during a crossing are geometrically consistent
    eij = np.sqrt(ecc_encounter[0]**2 + ecc_encounter[1]**2) #relative eccentricity 
 
    #calculates collosion probability, Pcol
    ln_lambda = 3
    lambdaa = (2*eij/e_esc)**2 * (1 + (eij**2/e_esc**2)) * (1/ln_lambda) #eq 12
    p_col = 1 - np.exp(-lambdaa) # eq 14, caps between [0,1]
 
    #use Monte Carlo approach to say whether or not a collision actually happens
    #Bernoulli sampling?
    merge_record = None #stores when a merge happens, stays None if nothing happens
    draw = np.random.uniform(0,1) #draws a random number to compare against p_col
    if draw < p_col: #merge event
        count = 0 #keeps track of rejection-sampling structure, and ensures while loop doesnt go forever
        while True: #keeps sampling until the eccentricity is consistent with the orbits actually overlapping 
            rayleigh_ecc = rayleigh(1 / np.sqrt(2), eij / e_esc)
            #mutate ecc[icross]/ecc[jcross] in place so each retry's max() compares against the
            #running (already-updated) eccentricity, matching the Fortran goto-loop behaviour
            ecc[icross] = max(rayleigh_ecc * np.sqrt(Mp[jcross]) / np.sqrt(Mp[icross] + Mp[jcross]) * e_esc, ecc[icross])
            ecc[jcross] = max(rayleigh_ecc * np.sqrt(Mp[icross]) / np.sqrt(Mp[icross] + Mp[jcross]) * e_esc, ecc[jcross])
 
            if np.sqrt(ecc_cross[0]**2 + ecc_cross[1]**2) / e_esc > 2.0 or count > 500:
                #unable to find a random draw that satisfies the condition, defaul to an orbital overlap of 0.1%
                #if orbits will never overlap OR took too many interations
                ecc[icross] = 1.001 * ecc_cross[0]
                ecc[jcross] = 1.001 * ecc_cross[1]
                break
            #check if epicycle amplitudes sum to at least aj-ai (they will overlap)
            #OVERLAP CONDITION
            elif ap[icross]*ecc[icross] + ap[jcross]*ecc[jcross] >= abs(ap[jcross] - ap[icross]):
                #yay it worked, eccentricities are already updated above
                break
            else: #if all else, try again with a new random number
                count += 1
 
        #velocity at collision, calculated from energy conservation
        v_c = collision_velocity(ap[icross:jcross+1], Mp[icross:jcross+1], Rp[icross:jcross+1], Ms, ecc[icross:jcross+1])

        #capture pre-merge state for the impact log, since merge_embryo overwrites Mp[icross:jcross+1] below
        #(target = larger/surviving body, impactor = smaller/destroyed body, as per merge_embryo convention)
        target_idx = icross if Mp[icross] >= Mp[jcross] else jcross
        impactor_idx = jcross if target_idx == icross else icross
        id_target, id_impactor = planet_id[target_idx], planet_id[impactor_idx] #carry index through all other parameters
        M_target_before, M_impactor_before = Mp[target_idx], Mp[impactor_idx]       
        
        
        #call merge_embryo function to update parameters for interacting pair
        #jcross+1 to include that planet in the interacting pair
        #print(f"[COLLISION] Planets {planet_id[icross]} and {planet_id[jcross]} merged")
        ap_merge, Mp_merge, ecc_merge, live_status_merge, atm_mass_fraction_merge, frac_lost = merge_embryo(ap[icross:jcross+1], Mp[icross:jcross+1], Rp[icross:jcross+1], Ms, ecc[icross:jcross+1], v_c, live_status[icross:jcross+1], impact_parameter, atm_mass_fraction[icross:jcross+1])
        #update system
        ap[icross:jcross+1] = ap_merge
        Mp[icross:jcross+1] = Mp_merge 
        ecc[icross:jcross+1] = ecc_merge 
        live_status[icross:jcross+1] = live_status_merge #smaller planet dies
        atm_mass_fraction[icross:jcross+1] = atm_mass_fraction_merge #should just apply to the target planet here

        #record impact velocity (v_c) + mass loss for each merger
        merge_record = {'id_target': id_target,'id_impactor': id_impactor,'M_target_before': M_target_before,
            'M_impactor_before': M_impactor_before,'M_merged_after': Mp[target_idx],'v_c': v_c,
            'atm_mass_loss_frac': frac_lost, 'a_final_AU': ap[target_idx] / 1.5e11,}

    else: #scattering event
        rayleigh_ecc = rayleigh(1.0 / np.sqrt(2.0), 0.0) #0 because no truncation at geometric overlap constraint as for merge case 
        #assuming energy equipartition
        #taken from fortran code, why is it not a max() anymore as before with merging?
        
        ecc[icross] = rayleigh_ecc * np.sqrt(Mp[jcross]) / np.sqrt(Mp[icross] + Mp[jcross]) * e_esc
        ecc[jcross] = rayleigh_ecc * np.sqrt(Mp[icross]) / np.sqrt(Mp[icross] + Mp[jcross]) * e_esc
 
        ilarge = icross if ecc[icross] >= ecc[jcross] else jcross #identify which planet was excited more
        ismall = jcross if ecc[icross] >= ecc[jcross] else icross
 
        if max(ecc[icross], ecc[jcross]) >= 1.0: #planet got bumped out
            #print(f"[EJECTION] Planet {planet_id[ismall]} was ejected")
            #planet with smaller excited eccentricity remains in the system, and orbital parameters are recalculated
            ap[ismall] = Mp[ismall] / (Mp[ismall] / ap[ismall] + Mp[ilarge] / ap[ilarge])
            #in K25 pt 2, this is ap[ismall] / aM - 1 
            ecc[ismall] = 1.0 - ap[ismall] / aM #mass-weighted mean, matches Fortran's aM (not the simple mean mma)
            live_status[ismall] = False
        else: #'normal' scattering conditions
            #'change in orbital separation is assumed to be equal to the sum of the excited epicycle amplitude'
            #db = delta_a essentially how much the orbit is shifted either in or out 
            #print(f"[SCATTERING] Planets {planet_id[icross]} and {planet_id[jcross]} scattered")
            db = ecc[icross] * ap[icross] + ecc[jcross] * ap[jcross] #delta b_ij eq 18
            ap[icross] = ap[icross] - Mp[jcross] / (Mp[icross] + Mp[jcross]) * db #eq 19
            ap[jcross] = ap[jcross] + Mp[icross] / (Mp[icross] + Mp[jcross]) * db #eq 20
 
    if ap[icross] < 0.0: #inner planet fell into star oops
        live_status[icross] = False #it's now dead
    
    return merge_record