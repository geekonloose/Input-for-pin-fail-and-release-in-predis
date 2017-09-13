'''
rf = pallet radious (meter)
dr = space between two grid point
alpha = 1/rho*Cp
input_data.dt_nu = spacing between two time points (seconds)
sp = no of grid point
time _step = no of time iteration
kc = clad conductance (W/m-k)
G = gap between clad and fuel (meter)
Q  = heat generation term
tc = clad thickness (meter)
T = temperature in fuel pin (K)

'''

import scipy, pylab, random
import numpy as np

np.set_printoptions(threshold=np.inf)
from matplotlib import pyplot as plt
import math
from inputs_for_sv_cracks import input_data
from temperature_profile import temp_profile
from crack_surface_area import crack
from u_sv_gb import ugr_area
from nu_modelling import nu_modelling

input_data()

T = temp_profile()

input_data.temp_T[0] = sum(T[1:11]) / 11
input_data.temp_T[1] = sum(T[12:22]) / 11
input_data.temp_T[2] = sum(T[23:33]) / 11

crack()

i1 = 0

for T in input_data.temp_T:
    [grain_area_single, D] = ugr_area(T)
    input_data.grainarea[i1] = grain_area_single
    input_data.Diff[i1] = D
    i1 = i1 + 1

'''
kr-89 = 3.75E-3
xe-137 = 3.01E-3
xe-138 = 8.17E-4
xe-135m = 7.55E-4
xe-134m = 7.35E-4
kr-87 = 1.515E-4
kr-83m = 1.05E-4
kr-88 = 6.73E-5
kr-85m = 4.415E-5
xe-135 = 8.12E-5
xe-133m = 3.6624E-6
xe-133 = 1.5424E-6
xe-131m = 6.7401E-7
kr-85 = 8.05E-9
'''

# print(input_data.isotopes)


# input_data.lmbda = np.array([1E-6, 1.53E-6, 9.26E-6, 8.12E-5, 8.91E-5, 4.30E-5, 6.78E-5, 8.37E-5,
# 1.52E-4, 8.2E-4, 7.56E-4, 8.18E-4, 3.02E-3])






# print('input_data.no_of_annuli',input_data.no_of_annuli)

for indx in range(input_data.lmbda.size):
    for anul in range(input_data.no_of_annuli):
        # print('s_v_c[%d]' % anul, input_data.s_v_c[anul])
        # print('input_data.grainarea[%d]' % anul, input_data.grainarea[anul])
        input_data.s_v_t[anul] = input_data.s_v_c[anul] + input_data.grainarea[anul]

        input_data.rba[anul] = input_data.s_v_t[anul] * (input_data.Diff[anul]  / input_data.lmbda[indx]) ** (0.5)
        # print('input_data.rba[%d]' % anul, input_data.rba[anul])
    input_data.rbf[indx] = sum(input_data.rba)

for indx in range(input_data.lmbda.size):
    if input_data.rbf[indx]>1:
        input_data.rbf[indx]=1

    print('input_data.rbf[%d]' % indx, input_data.rbf[indx])

'''
plt.plot(input_data.lmbda, input_data.rbf)
plt.semilogx()
plt.semilogy()
plt.show()
'''

# Program to find the release





'''
for i1 in range(input_data.lmbda.size):
    input_data.Na_inf[i1] = input_data.rbf[i1] * input_data.production[i1] / (input_data.lmbda[i1] + input_data.nu)
    input_data.Nb_inf[i1] = input_data.nu * input_data.Na_inf[i1] / (input_data.lmbda[i1] + input_data.sweep_to_volume)
    # print(input_data.isotopes[i1], (1 / (4E7)) * input_data.Nb_inf[i1])
'''
# for i1 in range(8,14):

# input_data.P[0] = 0
# input_data.rho = np.zeros([input_data.no_of_iteration])
# input_data.L = np.zeros([input_data.no_of_iteration])


temp_Na = np.zeros([input_data.lmbda.size, input_data.no_of_iteration])
temp_Nb = np.zeros([input_data.lmbda.size, input_data.no_of_iteration])
temp_P = np.zeros([input_data.lmbda.size, input_data.no_of_iteration])
# temp_P1 = np.zeros([input_data.lmbda.size,input_data.no_of_iteration])
temp_L = np.zeros([input_data.lmbda.size, input_data.no_of_iteration])
sum_of_Na = np.zeros([input_data.lmbda.size])
sum_of_Nb = np.zeros([input_data.lmbda.size])


for j1 in range(input_data.no_of_iteration - 1):
    if j1 < (input_data.pin_failure_time) and j1< input_data.reactor_shut_down_time:

        inv1 = (1 + input_data.lmbda * input_data.dt_nu) ** (-1)
        inv2 = (1 + input_data.lmbda * input_data.dt_nu + input_data.sweep_to_volume * input_data.dt_nu) ** (-1)

        input_data.Na = (input_data.rbf * input_data.dt_nu * input_data.production + input_data.Na) * inv1
        input_data.Nb = input_data.Nb * inv2

        temp_Na[8, j1] = input_data.Na[8]
        temp_Nb[8, j1] = input_data.Nb[8]

        # sum_of_Na = sum_of_Na + input_data.Na

        input_data.P1 = ((sum(input_data.Na) * 1.23E-22 * input_data.Tg) / (input_data.plenum_volume))  # 6.022E23 *

        temp_P[8, j1] = input_data.P1

    elif j1 >= input_data.pin_failure_time and j1<input_data.reactor_shut_down_time:
        print '2'
        input_data.constant = (input_data.rupture_dia**2 * input_data.area_of_rupture) / (
        64.0 * input_data.length_of_clad * input_data.mu)

        input_data.P1 = ((sum(input_data.Na) * 1.23E-22 * input_data.Tg) / (input_data.plenum_volume)) # 6.022E23 *


        if input_data.P1< input_data.P_ref:
            input_data.P1 = input_data.P_ref
        input_data.L = (input_data.constant * (input_data.P1 - input_data.P_ref))
        input_data.L = (1 / input_data.plenum_volume) * input_data.L
        print input_data.L
        inv1 = (1 + input_data.lmbda * input_data.dt_nu + input_data.L*input_data.dt_nu) ** (-1)
        inv2 = (1 + input_data.lmbda * input_data.dt_nu + input_data.sweep_to_volume * input_data.dt_nu) ** (-1)

        input_data.Na = (input_data.rbf * input_data.dt_nu * input_data.production + input_data.Na) * inv1
        input_data.Nb = (input_data.Nb + input_data.L * input_data.Na * input_data.dt_nu) * inv2



        temp_P[8, j1] = input_data.P1
        temp_L[8, j1] = input_data.L
        temp_Na[8, j1] = input_data.Na[8]
        temp_Nb[8, j1] = input_data.Nb[8]

    elif j1>=input_data.reactor_shut_down_time:

        print '3'
        input_data.constant = (input_data.rupture_dia ** 2 * input_data.area_of_rupture) / (
            64 * input_data.length_of_clad * input_data.mu)

        input_data.P1 = ((sum(input_data.Na) * 1.23E-22 * input_data.Tg) / (input_data.plenum_volume))  # 6.022E23 *

        if input_data.P1 < input_data.P_ref:
            input_data.P1 = input_data.P_ref
        input_data.L = (input_data.constant * (input_data.P1 - input_data.P_ref))
        input_data.L =  (1 / input_data.plenum_volume) * input_data.L

        inv1 = (1 + input_data.lmbda * input_data.dt_nu + input_data.L * input_data.dt_nu) ** (-1)
        inv2 = (1 + input_data.lmbda * input_data.dt_nu + input_data.sweep_to_volume * input_data.dt_nu) ** (-1)
        input_data.Na = (0.00+ input_data.Na) * inv1
        input_data.Nb = (input_data.Nb + input_data.L * input_data.Na * input_data.dt_nu) * inv2

        temp_P[8, j1] = input_data.P1
        temp_L[8, j1] = input_data.L
        temp_Na[8, j1] = input_data.Na[8]
        temp_Nb[8, j1] = input_data.Nb[8]
time_t = np.array([i4 for i4 in range(input_data.no_of_iteration)])
# print "L",input_data.L
plt.plot(input_data.dt_nu * time_t, temp_Na[8, :])
plt.show()

plt.plot(input_data.dt_nu * time_t, temp_Nb[8, :])
plt.show()

plt.plot(input_data.dt_nu * time_t, temp_P[8, :])
plt.show()

plt.plot(input_data.dt_nu * time_t, temp_L[8, :])
plt.show()


print sum(input_data.production)