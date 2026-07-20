import numpy as np
import matplotlib.pyplot as plt 
from astropy.table import Table
from astropy.io import ascii
import pdb
import toml
import os 
from constants import *
from scipy.stats import variation


plt.rcParams["font.size"] = 15
plt.rcParams['axes.linewidth'] = 2 
plt.rcParams["xtick.major.width"] = 1.5             
plt.rcParams["ytick.major.width"] = 1.5
plt.rcParams["xtick.minor.width"] = 1  
plt.rcParams["ytick.minor.width"] = 1  
plt.rcParams["xtick.major.size"] = 12  
plt.rcParams["ytick.major.size"] = 12  
plt.rcParams["xtick.minor.size"] = 7   
plt.rcParams["ytick.minor.size"] = 7  

with open('initialise.toml', 'r') as f:
    config = toml.load(f)

Ms = config['init_par']['Ms'] * M_sun 
ndisk = config.get('batch', {}).get('ndisk', 100) #number of systems run
directory = config['run_simulation']['save_directory']

#remaining planets in each system 
batch_summary = ascii.read(directory+'/batch_summary.csv', format = 'fixed_width')

def plot_tracks(directory):

    for n in range(ndisk):
        full_system = ascii.read(directory+f'/data/full_system_{n:02d}.csv', format = 'fixed_width')
        initial_N = N = config['init_par']['N'] #initial planets in each system

        for p in range(initial_N): 
            planet = full_system[full_system['id'] == p]
            if len(planet) == 0:
                continue
            t_yr = planet['t']/365/24/60/60
            a = planet['a_AU']
            e = planet['ecc']
    
            plt.grid(alpha = 0.5)
            line, = plt.plot(t_yr, a, label = f'Planet {p}')
            color = line.get_color()
    
            # Plot pericenter/apocenter boundaries and fill the orbital sweep region
            plt.plot(t_yr, a * (1.0 - e), linestyle=':', alpha=0.4, color=color)
            plt.plot(t_yr, a * (1.0 + e), linestyle=':', alpha=0.4, color=color)
            plt.fill_between(t_yr, a * (1.0 - e), a * (1.0 + e), alpha=0.1, color=color)

        plt.legend(ncol = 2, loc = 'upper right')
        plt.xlabel('Time (yr)')
        plt.ylabel('Radial range swept out by orbit (AU)')
        plt.xscale('log')
        plt.tight_layout()
        plt.savefig(directory+f'/figures/tracks/track{n:02d}.png', dpi = 300)
        plt.close()
 
plot_tracks(directory)

def plot_orbits(directory):

    #arrays to store information for all remaining planets in every system
    batch_planets = []
    batch_a = []
    batch_ecc = []
    batch_mp = []
    for n in range(ndisk):
        full_system = ascii.read(directory+f'/data/full_system_{n:02d}.csv', format = 'fixed_width')
        #pull remaining planets from each system
        batch = batch_summary[batch_summary['run_idx'] == n]
        n_survivors = int(batch['n_survivors'][0])
        #data for the planets left after simulation
        remaining_system = full_system[-n_survivors:]

        #plot orbits just for fun 
        fig, ax = plt.subplots(figsize=(7, 7))
        theta = np.linspace(0, 2 * np.pi, 300)
        c = ['palevioletred', 'steelblue', 'midnightblue', 'mediumseagreen', 'black']

        masses = remaining_system['Mp']/M_earth
        # scale marker sizes relative to masses in this system (area ~ mass, so use sqrt)
        size_min, size_max = 80, 600
        if masses.max() > masses.min():
            marker_sizes = size_min + (np.sqrt(masses) - np.sqrt(masses.min())) / \
                            (np.sqrt(masses.max()) - np.sqrt(masses.min())) * (size_max - size_min)
        else:
            marker_sizes = np.full(len(masses), (size_min + size_max) / 2)

        for a, e, mp, ms, f in zip(remaining_system['a_AU'], remaining_system['ecc'], masses, marker_sizes, range(len(remaining_system['ecc']))):
            r = a * (1 - e**2) / (1 + e * np.cos(theta))
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            ax.plot(x, y, lw=3, color=c[f], label=f'e = {round(e,3)}, Mp = {mp}' + r' $M_\oplus$')

            # mark planet position at periapsis (theta=0)
            r_peri = a * (1 - e)
            ax.scatter(r_peri, 0, s=ms, color=c[f], edgecolor='k', zorder=6)

        ax.plot(0, 0, marker='*', color='gold', markersize=15, markeredgecolor='orange', zorder=5)
        ax.set_aspect('equal')
        ax.set_xlabel('x (AU)')
        ax.set_ylabel('y (AU)')
        plt.grid()
        plt.legend(loc = 'upper right')
        plt.tight_layout()
        plt.savefig(directory+f'/figures/orbits/orbits_{n:02d}.png', dpi = 500)
        plt.close()


        batch_planets.append(len(remaining_system['id']))
        batch_a.append(remaining_system['a_AU'][0])
        batch_ecc.append(remaining_system['ecc'][0])
        batch_mp.append(remaining_system['Mp'][0]/M_earth)

plot_orbits(directory)

def mutual_hill(M1,M2,a1,a2,Ms):
    return ((M1 + M2)/(3 * Ms))**(1/3) * (a1 + a2)/2

def mean_a_ecc(M,a,e,Ms):
    M = np.asarray(M, dtype = float)
    a = np.asarray(a, dtype = float)
    e = np.asarray(e, dtype = float)

    Np = len(a)
    if Np < 2:
        raise ValueError('Need at least 2 remaining planets')
    
    order = np.argsort(a)
    M,a,e = M[order], a[order], e[order]

    b_H_terms = []
    e_H_terms = []

    for i in range(Np-1):
        r_H = mutual_hill(M[i], M[i+1], a[i], a[i+1], Ms)

        b_H_terms.append((a[i+1] - a[i]) / r_H)
        e_H_terms.append((e[i]*a[i] + e[i+1]*a[i+1]) / (2*r_H))

    b_H_mean = np.mean(b_H_terms)
    e_H_mean = np.mean(e_H_terms)

    return b_H_mean, e_H_mean

#initial center of mass for a singular system
def com_i(Mp,a):
    a_m = np.empty(len(a)) #for single system
    num = np.empty(len(a))
    den = np.empty(len(a))
    for i in range(len(a)):
        num[i] = Mp[i]*a[i]
        den[i] = Mp[i]

    a_M = sum(num)/sum(den)
    return a_M


def plot_stats(directory):
    # Arrays to accumulate information for all remaining planets across ALL systems (used in all_systems.png)
    batch_planets = []
    batch_a = []
    batch_ecc = []
    batch_mp = []
    
    # Pre-allocate arrays for per-system statistics (initialized to NaN to handle cases with < 2 survivors)
    std_a_run = np.full(ndisk, np.nan)
    std_M_run = np.full(ndisk, np.nan)
    com_run = np.full(ndisk, np.nan)
    b_H_run = np.full(ndisk, np.nan)
    e_H_run = np.full(ndisk, np.nan)
    
    for n in range(ndisk):
        full_system = ascii.read(directory+f'/data/full_system_{n:02d}.csv', format = 'fixed_width')
        # Pull remaining planets from each system
        batch = batch_summary[batch_summary['run_idx'] == n]
        n_survivors = int(batch['n_survivors'][0])
        # Data for the planets left after simulation
        remaining_system = full_system[-n_survivors:]
        # Extract system-specific attributes for THIS run to avoid cross-system accumulation
        run_a = np.asarray(remaining_system['a_AU'], dtype=float)
        run_ecc = np.asarray(remaining_system['ecc'], dtype=float)
        run_mp = np.asarray(remaining_system['Mp']/M_earth, dtype=float)
        run_mp_kg = np.asarray(remaining_system['Mp'], dtype=float)
        batch_planets.append(n_survivors)
        batch_a.extend(run_a)
        batch_ecc.extend(run_ecc)
        batch_mp.extend(run_mp)
        # Calculate statistics only for systems with sufficient survivors
        if n_survivors >= 1:
            com_run[n] = com_i(run_mp, run_a)
            if n_survivors >= 2:
                std_a_run[n] = variation(run_a, ddof = 0)
                std_M_run[n] = variation(run_mp, ddof = 0)
                b_H_run[n], e_H_run[n] = mean_a_ecc(run_mp_kg, run_a, run_ecc, Ms)


    fig,ax = plt.subplots(1,3, figsize = (15,5))
    ax[0].scatter(com_run, batch_planets, color = 'purple')
    ax[0].plot(np.nanmean(com_run), np.nanmean(batch_planets), marker='.', color='gold', markersize=15, markeredgecolor='black', zorder=5)
    ax[0].errorbar(np.nanmean(com_run), np.nanmean(batch_planets), yerr=np.nanstd(batch_planets), xerr=np.nanstd(com_run), ecolor = 'purple', capsize = 4, fmt='none')  
    ax[0].set_ylabel('Remaining planets')
    ax[0].set_xlabel(r'$\langle a_M \rangle$ (AU)')
    ax[0].grid(alpha = 0.5) 

    ax[1].scatter(e_H_run, b_H_run, color = 'seagreen')
    ax[1].plot(np.nanmean(e_H_run), np.nanmean(b_H_run), marker='.', color='midnightblue', markersize=15, markeredgecolor='black', zorder=5)
    ax[1].errorbar(np.nanmean(e_H_run), np.nanmean(b_H_run), yerr=np.nanstd(b_H_run), xerr=np.nanstd(e_H_run), ecolor = 'seagreen', capsize = 4, fmt='none') 
    ax[1].set_xlabel('Mean eccentricity (Hill-scaled)')
    ax[1].set_ylabel('Mean orbital separation (Hill-scaled)')
    ax[1].grid(alpha = 0.5)


    ax[2].scatter(std_a_run, std_M_run, color = 'coral')
    ax[2].plot(np.nanmean(std_a_run), np.nanmean(std_M_run), marker='.', color='gray', markersize=15, markeredgecolor='black', zorder=5)
    ax[2].errorbar(np.nanmean(std_a_run), np.nanmean(std_M_run), yerr=np.nanstd(std_M_run), xerr=np.nanstd(std_a_run), ecolor = 'coral', capsize = 4, fmt='none')
    ax[2].set_xlabel('Normalised a standard deviation (AU)')
    ax[2].set_ylabel('Normalised M standard deviation ($M_\oplus$)')
    ax[2].grid(alpha = 0.5)

    plt.tight_layout()
    plt.savefig(directory+'/figures/stats/all_stats.png', dpi = 300)
    plt.close()

    fig,ax = plt.subplots(1,3, figsize = (15,5))
    ax[0].scatter(batch_a, batch_ecc, color = 'midnightblue')
    ax[0].plot(np.mean(batch_a), np.mean(batch_ecc), marker='.', color='gold', markersize=15, markeredgecolor='black', zorder=5)
    ax[0].errorbar(np.mean(batch_a), np.mean(batch_ecc), np.std(batch_ecc), np.std(batch_a), ecolor = 'midnightblue', capsize = 4)
    ax[0].set_xlabel('a (AU)')
    ax[0].set_ylabel('ecc')
    ax[0].grid(alpha = 0.5)

    ax[1].scatter(batch_mp, batch_ecc, color = 'palevioletred')
    ax[1].plot(np.mean(batch_mp), np.mean(batch_ecc), marker='.', color='steelblue', markersize=15, markeredgecolor='black', zorder=5)
    ax[1].errorbar(np.mean(batch_mp), np.mean(batch_ecc), np.std(batch_ecc), np.std(batch_mp), ecolor = 'palevioletred', capsize = 4)
    ax[1].set_xlabel('M$_p$ ($M_\oplus$)')
    ax[1].grid(alpha = 0.5)

    ax[2].scatter(batch_a, batch_mp, color = 'forestgreen')
    ax[2].plot(np.mean(batch_a), np.mean(batch_mp), marker='.', color='purple', markersize=15, markeredgecolor='black', zorder=5)
    ax[2].errorbar(np.mean(batch_a), np.mean(batch_mp), np.std(batch_mp), np.std(batch_a), ecolor = 'forestgreen', capsize = 4)
    ax[2].set_xlabel('a')
    ax[2].set_ylabel('M$_p$ ($M_\oplus$)')
    ax[2].grid(alpha = 0.5)

    plt.tight_layout()
    plt.savefig(directory+'/figures/stats/all_systems.png', dpi = 500)
    plt.close()




plot_stats(directory)