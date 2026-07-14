import numpy as np
import matplotlib.pyplot as plt 
from astropy.table import Table
from astropy.io import ascii
import pdb
import toml
import os 

with open('initialise.toml', 'r') as f:
    config = toml.load(f)

directory = config['run_simulation']['save_directory']
os.makedirs(directory, exist_ok=True)

initial_system = ascii.read(directory+'/initial_system.csv', format = 'fixed_width')
final_system = ascii.read(directory+'/final_system.csv', format = 'fixed_width')
full_system = ascii.read(directory+'/full_system.csv', format = 'fixed_width')
initial_N = len(initial_system['a_AU'])
final_N = len(final_system['a_AU'])

fig,ax = plt.subplots(2,2)

ax[0,0].grid(alpha = 0.5)
ax[0,0].scatter(initial_system['a_AU'], range(initial_N), marker = '+', color = 'purple', label = 'Initial system')
ax[0,0].scatter(final_system['a_AU'], range(final_N), marker = 'x', color = 'green', label = 'Final system')
ax[0,0].set_xlabel('Semi-major axis (AU)')
ax[0,0].set_ylabel('Planets in system')

ax[0,1].grid(alpha = 0.5)
ax[0,1].scatter(initial_system['e'], range(initial_N), marker = '+', color = 'purple', label = 'Initial system')
ax[0,1].scatter(final_system['e'], range(final_N), marker = 'x', color = 'green', label = 'Final system')
ax[0,1].legend()
ax[0,1].set_xlabel('Eccentricity')

ax[1,0].grid(alpha = 0.5)
ax[1,0].scatter(initial_system['Mp']/5.9736e24, range(initial_N), marker = '+', color = 'purple', label = 'Initial system')
ax[1,0].scatter(final_system['Mp']/5.9736e24, range(final_N), marker = 'x', color = 'green', label = 'Final system')
ax[1,0].set_xlabel('Mass ($M_\oplus$)')
ax[1,0].set_ylabel('Planets in system')

ax[1,1].grid(alpha=0.5)
ax[1,1].scatter(initial_system['Rp']/6.371e6, range(initial_N), marker = '+', color = 'purple', label = 'Initial system')
ax[1,1].scatter(final_system['Rp']/6.371e6, range(final_N), marker = 'x', color = 'green', label = 'Final system')
ax[1,1].set_xlabel('Radius ($R_\oplus$)')

plt.tight_layout()
plt.savefig(directory+'/top_level.png', dpi = 300)
plt.close()


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
plt.savefig(directory+'/a_tracks.png', dpi = 300)
plt.close()
 


