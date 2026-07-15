import numpy as np
import matplotlib.pyplot as plt 
from astropy.table import Table
from astropy.io import ascii
import pdb
import toml
import os 

with open('initialise.toml', 'r') as f:
    config = toml.load(f)

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

def plot_stats(directory):

    #arrays to store information for all remaining planets in every system
    batch_planets = []
    batch_a = []
    batch_ecc = []
    for n in range(ndisk):
        full_system = ascii.read(directory+f'/data/full_system_{n:02d}.csv', format = 'fixed_width')
        #pull remaining planets from each system
        batch = batch_summary[batch_summary['run_idx'] == n]
        n_survivors = int(batch['n_survivors'])
        #data for the planets left after simulation
        remaining_system = full_system[-n_survivors:]

        #plot orbits just for fun 
        fig, ax = plt.subplots(figsize=(7, 7))
        theta = np.linspace(0, 2 * np.pi, 300)
        c = ['palevioletred', 'steelblue', 'midnightblue', 'mediumseagreen', 'black']
        for a, e, f in zip(remaining_system['a_AU'], remaining_system['ecc'], range(len(remaining_system['ecc']))):
            r = a * (1 - e**2) / (1 + e * np.cos(theta))
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            ax.plot(x, y, lw=3, color = c[f], label = f'e = {round(e,3)}')

        ax.plot(0, 0, marker='*', color='gold', markersize=15, markeredgecolor='orange', zorder=5)
        ax.set_aspect('equal')
        ax.set_xlabel('x (AU)')
        ax.set_ylabel('y (AU)')
        plt.grid()
        plt.legend()
        plt.tight_layout()
        plt.savefig(directory+f'/figures/orbits/orbits_{n:02d}.png', dpi = 500)
        plt.close()


        batch_planets.append(len(remaining_system['id']))
        batch_a.append(remaining_system['a_AU'][0])
        batch_ecc.append(remaining_system['ecc'][0])

    mean_N, std_N = np.mean(batch_planets), np.std(batch_planets)
    mean_a, std_a = np.mean(batch_a), np.std(batch_a)
    mean_ecc, std_ecc = np.mean(batch_ecc), np.std(batch_ecc)

    fig,ax = plt.subplots(1,2)

    ax[0].errorbar(mean_a, mean_N, std_N, std_a, fmt = 'o', color = 'steelblue', capsize = 3)
    ax[0].set_xlabel('Semi-major axis (AU)')
    ax[0].grid()
    ax[1].errorbar(mean_ecc, mean_N, std_N, std_ecc, fmt = 'o', color = 'palevioletred', capsize = 3)
    ax[1].set_xlabel('Eccentricity')
    ax[0].set_ylabel('Remaining planets')
    ax[1].grid()

    plt.tight_layout()
    plt.savefig(directory+'/figures/stats/all_stats.png', dpi = 500)
    plt.close()


plot_stats(directory)