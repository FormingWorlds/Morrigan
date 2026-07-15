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
os.makedirs(directory, exist_ok=True)
os.makedirs(directory + '/data', exist_ok = True)
os.makedirs(directory + '/figures', exist_ok = True)

def plot_tracks(directory):

    for n in range(ndisk):
        full_system = ascii.read(directory+f'/full_system_{n:02d}.csv', format = 'fixed_width')

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
        plt.savefig(directory+f'/figures/track{n:02d}.pdf', dpi = 300)
        plt.close()
 


