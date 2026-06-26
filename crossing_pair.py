import numpy as np 
from helper_functions.py import esc_ecc, hill_sphere_mutual 
from tau_cross.py import tau_vis, tau_col, interaction_wrapper 

def crossing_pair(ap, Mp, Rp, ecc, ecc_vec, g, beta, interact, N, t, t_ref): #identify crossing pair from triplet, return pair and t_event
    #evaluate tau_cross for all planets
    Tcross = np.full(N - 1, 1e20) #initial array to eventually store the predicted crossing times for every pair, initialised at 'never'
    Naffect = np.ones(N, dtype=int) #the interacting planets, used to compute K factor in eq 5
    interacting_pair = np.arange(N - 1) #stores which pair of the triplet interact

    #calculate current orbital separation
    bKmin = 2.5 #result from N-body simulation, Kokubo et al 2025
    bKij = np.empty(N - 1)
    for i in range(N-1):
        a_mean = 0.5 * (ap[i] + ap[i+1])
        #maximum separation before ejection
        rKij = esc_ecc(a_mean, Mp[i], Mp[i+1], Rp[i], Rp[i+1]) * a_mean #physical separation
        bKij[i]  = (ap[i+1]-ap[i])/rKij #current gap relative to physical separation

    #2-body case
    if N == 2:
        aM = (Mp[0]*ap[0] + Mp[1]*ap[1]) / (Mp[0] + Mp[1]) #mean semi-major axis
        h = hill_sphere_mutual(Mp[0]+Mp[1], aM) / aM
        #checking for stability again here
        EJbef = 5.0/8.0*(ecc[0]**2 + ecc[1]**2)/h**2 - 3.0/8.0 * ((ap[0]-ap[1])/(h*aM))**2 + 4.5
        if EJbef > 0.0: #system is stable
            return 0, 1.5 * t #no crossing, punts event forward indefinitely to end simulation

    #calculate orbit crossing time for each pair with the original eccentricities
    for i in range(N - 2):
        if not all(interact[i:i+3]):
            continue #skip triplets where a planet is not live or interacting 
        
        #determine 'packed planets'
        group = np.zeros(N, dtype=bool)
        #true if conditions are satisfied, otherwise false
        group[:i+1] = (bKij[:i+1] < bKmin) & interact[:i+1]
        group[i+1:] = (bKij[i:N-1] < bKmin) & interact[i+1:]

        #'let up to 2 neighboring planets inside and outside the triplet affect the number of resonances'
        #find indices of planets to the left (i_in) and right (i_out) of the triplet
        non_group_in = np.where(~group[:i+1])[0] #indices if planets not in the triplet, flips boolean
        if len(non_group_in) == 0 or all(group[:i+1]): #if there are no planets, or all of them are packed
            i_in = 0
        else: #otherwise, i_in is the index of the planet closest to inner edge of triplet
            i_in = min(non_group_in[-1] + 1, i)
        
        #similarly for outer planet
        non_group_out = np.where(~group[i+1:])[0]
        if len(non_group_out) == 0 or all(group[i+1:]):
            i_out = Np - 1 #if everything is packed, the group extends to the end of the system
        else: #index of planet closest to outer edge of the triplet
            i_out = max(non_group_out[0] + i, i + 1)

        #numerical factor to account for the assumption that resonance density in a N-body system is K times larger than in the 3-planet case
        #scales up 3-planet system to N-planet system 
        Naffect_val = max(i_out - i_in + 1, 3) #eq 5
        Naffect[i] = Naffect_val

        #chooses interacting pair from the triplet 
        d1 = (1.0 - ecc[i+1])*ap[i+1] - (1.0 + ecc[i])*ap[i] #eq 4, closest physical distance so not normalised 
        d2 = (1.0 - ecc[i+2])*ap[i+2] - (1.0 + ecc[i+1])*ap[i+1]
        pair_index = i if d1 <= d2 else i+1
        interacting_pair[i] = pair_index
        
        #eccentricity at overlap eq6
        ecc_cross = [np.sqrt(Mp[pair_index+1]) * abs(ap[pair_index+1]-ap[pair_index]) / (np.sqrt(Mp[pair_index+1])*ap[pair_index] + np.sqrt(Mp[pair_index])*ap[pair_index+1]),np.sqrt(Mp[pair_index]) * abs(ap[pair_index+1]-ap[pair_index]) / (np.sqrt(Mp[pair_index+1])*ap[pair_index] + np.sqrt(Mp[pair_index])*ap[pair_index+1])]
        
        #as before, the eccentricity can't be lower than eccentricity actually required for an overlap
        ecc_cross[0] = max(ecc[pair_index], ecc_cross[0]) #eq 23
        ecc_cross[1] = max(ecc[pair_index+1], ecc_cross[1])
        
        #calculate interaction timescales for every triplet, calling timescale functions from above
        viscous_timescale = tau_vis(ap[pair_index:pair_index+2], Mp[pair_index:pair_index+2], Rp[pair_index:pair_index+2], ecc_cross)
        collision_timescale = tau_col(ap[pair_index:pair_index+2], Mp[pair_index:pair_index+2], Rp[pair_index:pair_index+2], ecc_cross)
        
        #correction with secular perturbations, eq 21
        ecc_dmy = [np.sqrt(np.sum(ecc_vec[i, :]**2)),np.sqrt(np.sum(ecc_vec[i+1, :]**2)),np.sqrt(np.sum(ecc_vec[i+2, :]**2))]
        
        #passes the TRIPLETS ap, Mp, ecc to the wrapper function to calculate crossing timescale 
        crossing_timescale = interaction_wrapper(ap[i:i+3], Mp[i:i+3], ecc_dmy, Naffect_val)
        #predicted crossing timescale is 'current' time t + predicted next interaction + duration of the interaction
        Tcross[i] = t + crossing_timescale + min(viscous_timescale, collision_timescale)
        
    idmy_min = np.argmin(Tcross) #index of minimum crossing time
    t_event = Tcross[idmy_min]
    icross = interacting_pair[idmy_min] 

    return icross, t_event #inner planet index, event time