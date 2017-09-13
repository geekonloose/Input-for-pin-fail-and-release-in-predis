'''
rf = pallet radious (meter)
dr = space between two grid point
alpha = 1/rho*Cp
dt = spacing between two time points (seconds)
sp = no of grid point
time _step = no of time iteration
kc = clad conductance (W/m-k)
G = gap between clad and fuel (meter)
Q  = heat generation term
tc = clad thickness (meter)
T = temperature in fuel pin (K)

'''



import numpy as np
from matplotlib import pyplot as plt
import math

#%%input constants
rf = (5.52/2)*1E-3
dr = 1E-4  # mesh spacing
alphat = 1E-6  # alpha only contains 1/rho*Cp
dt_max = (dr ** 2) / (2 * alphat)  # stability condition
dt = dt_max / 10  # time step size
sp = int(rf/dr) # no of spatial points
time_step = int(1e4)  # time points
kc = 29.0  # clad conductance
kg = 58.22E-3  # conductance in gap between clad and gap
G = 427E-6  # gap
hg = kg / G  # gap conductance
linear_heat_rate = 56E3
Qe = linear_heat_rate/(3.14*rf**2)

#%% Parameter for calculation of k
k = [2 for i in range(sp + 1)]  # initial array of k
x = 2.0 - 1.97  # 1.97 is O/M ratio
A = 2.85 * x + 0.035
B = -0.715 * x + 0.286
space = [i for i in range(1, sp)]

def main():
    T  = temp_profile()
    plt.plot(  space, T[1:  sp],color = 'blue')
    plt.xlabel('radial distance(mm)',color = 'blue',size='14')
    plt.ylabel('Temperature(K)',color = 'blue',size='14')
    plt.show()
    plt.savefig('temperature_profile_contact1')








def temp_profile():


    T = [397 for i in range(  sp + 1)]

    # fuel surface temperature
    T[-1] = 883 + 273.0

    for j in range(0,   time_step):
        for i in range(1,   sp):
            T[1] = T[2]
            t = T[i] / 1000.0
            k[i] = 1.158*((1.0/float(  A+  B*  t)) + (6400.0/(float(  t)**(2.5))* math.exp(-16.35/float(  t))))
            T[i] = (T[i] +((  alphat *   dt) / (i * (  dr ** 2)))
* (((2 * i + 1) / 2) *   k[i+1] * (T[i + 1] - T[i]) - ((2 * i - 1) / 2) *   k[i-1] * (T[i] - T[i - 1])) + (  Qe *   alphat *   dt))



    # Tc1 = T[-1] - (Qe * rf) / (hg * 2)
    # Tc2 = Tc1 - ((Qe * rf * tc) / (kc * 2))
    # pallet_temp = [T[-1, -1], Tc1, Tc2]
    # #print('pallet temp ', pallet_temp)

    return (T)

if __name__=="__main__":
    main()