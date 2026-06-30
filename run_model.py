import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
import pdb 
from astropy.table import Table
from astropy.io import ascii
import os 
import time 

#import functions
from helper_functions import * 
from tau_cross import * 
from merge_embryo import * 
from secular_solution import * 
from crossing_pair import * 
from orbit_cross_K25 import * 
from sort_planet import *

#####CONSTANTS######
G = 6.67e-11 #m^3kg^-1s^-2
M_sun = 1.9892e30 #kg
M_earth = 5.9736e24 #kg

####VARIABLES TO RUN SIMULATION####
t = 0.0
t_ref = 0.0
t_event = 0.0
flag_event = 1 
a_min = 0.005 * 1.5e11 #defines when a planet has fallen into the star 
max_time = 5e9*365*24*60*60 #evolution time (seconds)
save_directory = 'test'
os.makedirs(save_directory, exist_ok=True) #creates directory if it doesnt already exist
####INITIAL PARAMETERS######
N = 20 #number of planets
e = 0.01
Mp = 0.5*M_earth #planet mass (relative to Mearth)
Ms = 1*M_sun #stellar mass (relative to Msun)
rho_p = 5500 #planet density kg/m^3
#np.random.seed(1) #for reproducability

####ALLOCATE PARAMETERS FOR THE SYSTEM###

def allocate_a(M):
    a = np.empty(N)
    a[0] = 0.1 #AU
    for i in range(1,N): #allocate initial semi-major axes
        a_previous = a[i-1] #starting semi-major axis
        hill = hill_sphere(a_previous,M)
        a[i] = a_previous + 10*hill #planets are spaced out by 10 hill radii
    return a*1.5e11 #convert to [m] to stay in SI!

#actually initialising system here with arrays for every parameter
a = allocate_a(Mp)
ecc = np.full(N,e)
densities = np.full(N, rho_p)
masses = np.full(N, Mp)
live_status = np.ones(N, dtype = bool) #set initial status of planets, all are live by definition at the start
interact = np.ones(N, dtype = bool) #stores the indices of which planets are participating in an event
Rp = np.array([planet_radius(i, j) for i,j in zip(masses,densities)])
planet_id = np.arange(N) #persistent id for a particular planet to track its evolution

parameter_names = ['id','a_AU','e','Mp','Rp','live_status']
system_information = Table([planet_id,a/1.5e11,ecc,masses,Rp,live_status],names = parameter_names)
ascii.write(system_information, save_directory+'/initial_system.csv', format = 'fixed_width', overwrite = True) 
#store initial system information

#timestep when not at or during an event
def time_step(t, t_event):
    dt = 0.1*(t+1.0e2*365*24*60*60) #timestep in seconds
    dt = min(dt,abs(t_event - t) + 1.0)
    return dt

#for tracking individual planet trajectories
planet_id = np.arange(len(masses))  # persistent IDs, assigned once at t=0
next_id = len(masses)

history = []
#stores timestep information about the system
def snapshot(t, a, masses, ecc, Rp, live_status, planet_id, N, event=False):
    history.append({'t': t, 'id': planet_id[:N].copy(), 'a': a[:N].copy()/1.5e11, 'masses': masses[:N].copy(),
        'ecc': ecc[:N].copy(),'Rp': Rp[:N].copy(),'live_status': live_status[:N].copy(),'event': event,})


output_interval = max_time / 1000.0  #when not at an event, store information every 1000 step
next_output = 0.0

# initial snapshot at t=0
snapshot(t, a, masses, ecc, Rp, live_status, planet_id, N, event=False)
next_output += output_interval

start = time.time()
#run simulation
while t <= max_time and N > 1:
    if flag_event == 1: #only recompute secular solution and crossing pair when something has changed
        a, masses, ecc, Rp, live_status, interact, densities, planet_id = sort_planet(a, masses, ecc, Rp, live_status, interact, densities, planet_id)
        N = len(a) #number of planets changes after an event!
        ecc_vec, g, beta = secular_solution(a, masses, ecc, Rp, N)
        t_ref = t #time for crossing_pair

        icross, t_event = crossing_pair(a, masses, Rp, ecc, ecc_vec, g, beta, interact, N, t, t_ref)
        flag_event = 0 #event done, do not recalculate secular/crossing otherwise 

    dt = time_step(t, t_event)
    t += dt #adjust time to account for event
    
    #propagate secular eccentricities
    h_t = np.zeros(N)
    k_t = np.zeros(N)
    for i in range(N):
        h_t[i] = np.sum(ecc_vec[i, :] * np.sin(g * (t - t_ref) + beta)) # eq A10
        k_t[i] = np.sum(ecc_vec[i, :] * np.cos(g * (t - t_ref) + beta)) 
    ecc = np.sqrt(h_t**2 + k_t**2) # eq 3 finally!
    
    #check for crossings/close encounters
    if t >= t_event:
        flag_event = 1
        orbit_cross_K25(a, masses, Rp, ecc, interact, live_status, N, planet_id, icross)
        snapshot(t, a, masses, ecc, Rp, live_status, planet_id, N, event=True) #capture state of system right after event
    
    #update planet radius
    Rp = np.array([planet_radius(masses[i], densities[i]) for i in range(N)])
    
    #remove planets too close to the star
    for i in range(N):
        if (1.0 - ecc[i]) * a[i] < a_min:
            live_status[i] = False
            flag_event = 1

    if t >= next_output:
        snapshot(t, a, masses, ecc, Rp, live_status, planet_id, N, event=False)
        next_output += output_interval

snapshot(t, a, masses, ecc, Rp, live_status, planet_id, N, event=False) #final system information

def data_to_table(history):
    t_col, id_col, a_col, m_col, e_col, rp_col, alive_col, event_col = [], [], [], [], [], [], [], []
    for h in history:
        n = len(h['id'])
        t_col.extend([h['t']] * n)
        id_col.extend(h['id'])
        a_col.extend(h['a'])
        m_col.extend(h['masses'])
        e_col.extend(h['ecc'])
        rp_col.extend(h['Rp'])
        alive_col.extend(h['live_status'])
        event_col.extend([h['event']] * n)
    return Table([t_col, id_col, a_col, m_col, e_col, rp_col, alive_col, event_col],names=['t', 'id', 'a_AU', 'Mp', 'ecc', 'Rp', 'live_status', 'event'])

ascii.write(data_to_table(history), save_directory+'/full_system.csv', format = 'fixed_width', overwrite = True)

#save final system information
system_information = Table([planet_id, a/1.5e11, ecc, masses, Rp, live_status], names=parameter_names)
ascii.write(system_information, save_directory+'/final_system.csv', format='fixed_width', overwrite=True)

end = time.time()
print(f'Runtime = ',round((end-start)/60, 3))