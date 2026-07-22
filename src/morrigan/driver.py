"""
!!! info "`driver.py`"
    Runs a dynamical evolution model for [ndisk] systems with [N] planets each and saves data in .csv files
    Author(s): Anna Grace Ulses
"""

import argparse
import os
import time
from functools import partial
from multiprocessing import Pool, cpu_count

import numpy as np
import toml
from astropy.io import ascii
from astropy.table import Table

#import functions and constants
from morrigan.constants import M_earth, M_sun, au2m, gyr2sec
from morrigan.crossing_pair import crossing_pair
from morrigan.helper_functions import hill_sphere, planet_radius
from morrigan.orbit_cross_K25 import orbit_cross_K25
from morrigan.secular_solution import secular_solution
from morrigan.sort_planet import sort_planet

#name of the settings file read when none is given on the command line
DEFAULT_CONFIG = 'initialise.toml'


def read_config(config_path=DEFAULT_CONFIG):
    """
    Read a settings file

    Parameters
    ----------
    config_path : str
        Path to the .toml settings file

    Returns
    -------
    config : dict
        Parsed settings
    """
    with open(config_path, 'r') as f:
        return toml.load(f)


####ALLOCATE PARAMETERS FOR THE SYSTEM###

def allocate_a(N,Ms,masses,inner_edge):
    a = np.empty(N)
    a[0] = inner_edge #AU
    for i in range(1,N): #allocate initial semi-major axes
        a_previous = a[i-1] #starting semi-major axis
        hill = hill_sphere(a_previous,masses[i-1]+masses[i],Ms) #mutual hill radius of the adjacent pair
        a[i] = a_previous + 10*hill #planets are spaced out by 10 hill radii
    return a*au2m #convert to [m] to stay in SI!

#timestep when not at or during an event
def time_step(t, t_event):
    dt = 0.1*(t+1.0e2*365*24*60*60) #timestep in seconds
    dt = min(dt,abs(t_event - t) + 1.0)
    return dt

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

def run_once(run_idx, config):

    #import settings from .toml file
    base_seed = config['run_simulation'].get('random_seed', 0)
    np.random.seed(base_seed + run_idx) #unique seed for each planetary system to reproduce individual results exactly

    t = config['run_simulation']['t']
    t_ref = config['run_simulation']['t_ref']
    t_event = config['run_simulation']['t_event']
    flag_event = config['run_simulation']['flag_event'] 

    a_min = config['run_simulation']['a_min'] * au2m #defines when a planet has fallen into the star 
    max_time = config['run_simulation']['max_time'] * gyr2sec #evolution time gyr converted to seconds
    save_directory = config['run_simulation']['save_directory']
    os.makedirs(save_directory, exist_ok=True) #creates directory if it doesnt already exist
    #sub-directories for data tables and figures respectively
    os.makedirs(save_directory+'/data', exist_ok = True)
    os.makedirs(save_directory+'/data/mergers', exist_ok = True)
    os.makedirs(save_directory+'/data/full_systems', exist_ok = True)
    os.makedirs(save_directory+'/data/survivors', exist_ok = True)
    os.makedirs(save_directory+'/figures', exist_ok = True)
    os.makedirs(save_directory + '/figures/tracks', exist_ok = True)
    os.makedirs(save_directory + '/figures/stats', exist_ok = True)
    os.makedirs(save_directory + '/figures/orbits', exist_ok = True)

    N = config['init_par']['N'] #number of planets
    e = config['init_par']['e'] #initial eccentricity
    impact_angle = config['init_par']['impact_angle']
    impact_parameter = np.sin(np.deg2rad(impact_angle)) #impact parameter =  sin(impact_angle)
    masses = np.array(config['init_par']['Mp']) * M_earth
    atm_mass_fraction = np.array(config['init_par']['atm_mass_fraction']) #mass fraction, solid mass implicitly is Mp - (Mp * atm_mass_fraction)

    if N != len(masses) or N != len(atm_mass_fraction):
        print('Number of planets = ', N)
        print('Masses allocated = ', len(masses))
        print('Atmospheres allocated = ', len(atm_mass_fraction))
        raise ValueError('Initial masses and/or atmosphere mass fractions and number of planets in system are mismatched!')

    Ms = config['init_par']['Ms'] * M_sun #stellar mass (relative to Msun)
    rho_p = config['init_par']['rho_p'] #planet density kg/m^3
    inner_edge = config['init_par']['inner_edge'] #orbit of the innermost planet (AU)

    #actually initialising system here with arrays for every parameter
    a = allocate_a(N,Ms,masses,inner_edge)
    ecc = np.full(N,e)
    densities = np.full(N, rho_p)
    live_status = np.ones(N, dtype = bool) #set initial status of planets, all are live by definition at the start
    interact = np.ones(N, dtype = bool) #stores the indices of which planets are participating in an event
    Rp = np.array([planet_radius(i, j) for i,j in zip(masses,densities)])
    planet_id = np.arange(N) #persistent id for a particular planet to track its evolution and what events it participates in
                                            
    history = []
    mergers = [] #specifically stores info about merge events, one row is one merge
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
            a, masses, ecc, Rp, atm_mass_fraction, live_status, interact, densities, planet_id = sort_planet(a, masses, ecc, Rp, atm_mass_fraction, live_status, interact, densities, planet_id)
            N = len(a) #number of planets changes after an event!
            if N <= 1:
                break 
            ecc_vec, g, beta = secular_solution(a, masses, ecc, Rp, Ms, N)
            t_ref = t #time for crossing_pair
            #identify indices of planetary pair that cross, and time of crossing (event)
            icross, t_event = crossing_pair(a, masses, Rp, Ms, ecc, ecc_vec, g, beta, interact, N, t, t_ref)
            flag_event = 0 #event done, do not recalculate secular/crossing otherwise 

        dt = time_step(t, t_event)
        t += dt #adjust time to account for event duration
    
        #propagate secular (long-term) eccentricities
        h_t = np.zeros(N)
        k_t = np.zeros(N)
        for i in range(N):
            h_t[i] = np.sum(ecc_vec[i, :] * np.sin(g * (t - t_ref) + beta)) # eq A10
            k_t[i] = np.sum(ecc_vec[i, :] * np.cos(g * (t - t_ref) + beta)) 
        ecc = np.sqrt(h_t**2 + k_t**2) # eq 3 finally!
    
        #check for next crossings/close encounters
        if t >= t_event:
            flag_event = 1
            merge_record = orbit_cross_K25(a, masses, Rp, Ms, atm_mass_fraction, impact_parameter, ecc, interact, live_status, N, planet_id, icross)
            if merge_record is not None: #None for scattering/ejection events, ONLY for mergers
                merge_record['t'] = t
                mergers.append(merge_record)
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

    snapshot(t, a, masses, ecc, Rp, live_status, planet_id, N, event=False) #final system snapshot

    ascii.write(data_to_table(history), os.path.join(save_directory+'/data/full_systems', f'full_system_{run_idx:02d}.csv'), format = 'fixed_width', overwrite = True)
    #write out impact velocities + resultant atmospheric mass loss for every merger in this run
    merger_cols = ['t', 'id_target', 'id_impactor', 'M_target_before', 'M_impactor_before',
                   'M_merged_after', 'v_c', 'atm_mass_loss_frac', 'a_final_AU']
    if mergers: #at least one merger happened this run (unlikely that there are no mergers)
        merger_table = Table(rows=mergers, names=merger_cols)
    else: #keep the file schema consistent even for runs with zero mergers
        merger_table = Table(names=merger_cols, dtype=[float]*len(merger_cols))
    ascii.write(merger_table, os.path.join(save_directory+'/data/mergers', f'mergers_{run_idx:02d}.csv'), format = 'fixed_width', overwrite = True)

    #save the final surviving planets for this run: id, mass, semi-major axis, eccentricity, remaining atmosphere fraction
    survivor_mask = live_status.astype(bool)
    survivors_table = Table([planet_id[survivor_mask], masses[survivor_mask], a[survivor_mask] / au2m, ecc[survivor_mask], atm_mass_fraction[survivor_mask]],
        names=['id', 'Mp', 'a_AU', 'ecc', 'atm_mass_fraction'])
    #atm_mass_fraction here is the fraction of the planet's mass that is atmosphere after all mergers
    ascii.write(survivors_table, os.path.join(save_directory+'/data/survivors', f'survivors_{run_idx:02d}.csv'), format = 'fixed_width', overwrite = True)

    end = time.time()
    runtime = round((end-start), 3)
    return {'run_idx': run_idx, 'runtime_s': runtime, 'n_survivors': int(np.sum(live_status))}

def main(config_path):
    """
    Run every system described by a settings file

    Parameters
    ----------
    config_path : str
        Path to the .toml settings file
    """
    #each disk is initialised with the same conditions
    config = read_config(config_path)

    #number of systems to run (defaults to a single run, unless specified)
    ndisk = config.get('batch', {}).get('ndisk', 1)

    os.makedirs(config['run_simulation']['save_directory'], exist_ok=True)

    if ndisk <= 1:
        #single run so skip multiprocessing entirely, allows for debugging
        start = time.time()
        result = run_once(0, config)
        end = time.time()
 
        summary = Table(rows=[result], names=['run_idx', 'runtime_s', 'n_survivors'])
        ascii.write(summary, os.path.join(config['run_simulation']['save_directory'], 'batch_summary.csv'),
                    format='fixed_width', overwrite=True)
        print(f'Ran 1 system in {round(end - start, 3)}s')
 
    else:
        #number of cores to use, leaves one free by default
        nproc = config.get('batch', {}).get('nproc', max(1, cpu_count() - 1))
 
        worker = partial(run_once, config=config)
 
        start = time.time()
        with Pool(processes=nproc) as pool:
            results = pool.map(worker, range(ndisk))
        end = time.time()
 
        #high-level statistics for each system (remaining planets, ids, etc)
        summary = Table(rows=results, names=['run_idx', 'runtime_s', 'n_survivors'])
        ascii.write(summary, os.path.join(config['run_simulation']['save_directory'], 'batch_summary.csv'),
                    format='fixed_width', overwrite=True)
        print(f'Ran {ndisk} systems in {round(end - start, 3)}s')


def cli():
    """
    Read the settings-file path from the command line, then run

    Kept separate from main() so that importing Morrigan and calling
    main() from another program never inspects that program's own
    command line.
    """
    parser = argparse.ArgumentParser(description='Run the Morrigan giant-impact model')
    parser.add_argument('-c', '--config', default=DEFAULT_CONFIG,
                        help='path to the .toml settings file')
    main(parser.parse_args().config)


if __name__ == '__main__':
    cli()

