import numpy as np 

#####CONSTANTS######
G = 6.67e-11 #m^3kg^-1s^-2
M_sun = 1.9892e30 #kg
M_earth = 5.9736e24 #kg

Ms = 1*M_sun #stellar mass (relative to Msun)

def kepler_period(Mp,a): #period of planetary orbit, used to calculate tau_cross
    P_squared = (4*np.pi**2*a**3)/(G*(Mp+Ms))
    return np.sqrt(P_squared)

def rayleigh(sigma, xmin):
    Umin = 1.0 - np.exp(-0.5 * xmin**2 / sigma**2)
    dum = np.random.uniform(Umin, 1.0 - 1e-10)
    return sigma * np.sqrt(-2.0 * np.log(1.0 - dum))

def esc_ecc(M1,M2,R1,R2,a):
    num = np.sqrt(2*G*(M1+M2)/(R1+R2))
    denom = np.sqrt((G*Ms)/a)
    return num/denom

def planet_radius(mass,density):
    return ((3*mass)/(4 * np.pi*density))**(1/3)

def hill_sphere(a_i,M):
    return a_i * ((M) / (3 * Ms))**(1/3) #mutual hill radius for adjacent planets