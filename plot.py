import numpy as np
import matplotlib.pyplot as plt 
from astropy.table import Table
from astropy.io import ascii

directory = 'jun29'
initial_system = ascii.read(directory+'/initial_system.csv', format = 'fixed_width')
final_system = ascii.read(directory+'/final_system.csv', format = 'fixed_width')
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
ax[1,1].scatter(initial_system['Rp']/6371e3, range(initial_N), marker = '+', color = 'purple', label = 'Initial system')
ax[1,1].scatter(final_system['Rp']/6371e3, range(final_N), marker = 'x', color = 'green', label = 'Final system')
ax[1,1].set_xlabel('Radius ($R_\oplus$)')


plt.tight_layout()
plt.savefig(directory+'/test.png', dpi = 300)
