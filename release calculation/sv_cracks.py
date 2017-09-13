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
import csv
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import math
import scipy, pylab, random
import numpy as np

np.set_printoptions(threshold=np.inf)

import math
from inputs_for_sv_cracks import input_data
from temperature_profile import temp_profile
from crack_surface_area import crack
from u_sv_gb import ugr_area
from nu_modelling import nu_modelling

# input_data()


T = temp_profile()



input_data.temp_T[0] = sum(T[2:7]) / 5
input_data.temp_T[1] = sum(T[7:12]) / 5
input_data.temp_T[2] = sum(T[12:17]) / 5
input_data.temp_T[3] = sum(T[17:22]) / 5
input_data.temp_T[4] = sum(T[22:28]) /5


crack()
print (input_data.temp_T)

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
xe-135 = 2.12E-5
xe-133m = 3.6624E-6
xe-133 = 1.5424E-6
xe-131m = 6.7401E-7
kr-85 = 2.05E-9
'''

# print(input_data.isotopes)


# input_data.lmbda = np.array([1E-6, 1.53E-6, 9.26E-6, 2.12E-5, 2.91E-5, 4.30E-5, 6.78E-5, 8.37E-5,
# 1.52E-4, 2.2E-4, 7.56E-4, 8.18E-4, 3.02E-3])






# print('input_data.no_of_annuli',input_data.no_of_annuli)

for indx in range(input_data.lmbda.size):
    for anul in range(input_data.no_of_annuli):
        # print('s_v_c[%d]' % anul, input_data.s_v_c[anul])
        # print('input_data.grainarea[%d]' % anul, input_data.grainarea[anul])
        input_data.s_v_t[anul] = input_data.s_v_c[anul] + input_data.grainarea[anul]

        input_data.rba[anul] = input_data.s_v_t[anul] * (input_data.Diff[anul]  / input_data.lmbda[indx]) ** (0.5)
        # print('input_data.rba[%d]' % anul, input_data.rba[anul])
    input_data.rbf[indx] = sum(input_data.rba)
    print('{}'.format(input_data.isotopes[indx]), input_data.rbf[indx])

plt.plot(input_data.lmbda, input_data.rbf,color='green',label = 'calculated')
plt.legend()
plt.semilogx()
plt.semilogy()
plt.show()

# %% Origin data

x = []
y = []
with open('RB data.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    for row in readCSV:
        try:
            x.append(row[0])
            y.append(row[1])
        except:
            pass

x1 = x[1:-1]
x1.append(x[-1])


x1 = list(map(float, x1))
y = list(map(float,y))
x1 = np.array(x1)
y = np.array(y)

def func(x1,a,b):
    return a* x1** b


params = curve_fit(func, x1, y)
[a, b] = params[0]
plt.plot(x1,y,color='red')
plt.semilogx()
plt.semilogy()
y = a* x1** b
plt.plot(x1,y,color='blue',ls='--',label = 'calculated by koo')
plt.legend()
plt.xlim(1E-7,0.01)
plt.ylim(1E-6,1)
plt.xlabel('lambda',color = 'blue',size = '13')
plt.ylabel('R/B',color = 'blue',size = '13')
plt.semilogx()
plt.semilogy()
plt.show()
plt.savefig('result')






'''
# Program to find the release






# for i1 in range(input_data.lmbda.size):
#     input_data.Na_inf[i1] = input_data.rbf[i1] * input_data.production[i1] / (input_data.lmbda[i1] + input_data.nu)
#     input_data.Nb_inf[i1] = input_data.nu * input_data.Na_inf[i1] / (input_data.lmbda[i1] + input_data.sweep_to_volume)
# print(input_data.isotopes[i1], (1 / (4E7)) * input_data.Nb_inf[i1])





# for i1 in range(13,14):

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

input_data.constant = (input_data.rupture_dia ** 2 * input_data.area_of_rupture) / (
    64 * input_data.length_of_clad * input_data.mu)

for j1 in range(input_data.no_of_iteration - 1):
    if (input_data.P1<2E6):

        inv1 = (1 + input_data.lmbda * input_data.dt_nu) ** (-1)
        inv2 = (1 + input_data.lmbda * input_data.dt_nu + input_data.sweep_to_volume * input_data.dt_nu) ** (-1)

        input_data.Na = (input_data.rbf * input_data.dt_nu * input_data.production + input_data.Na) * inv1
        input_data.Nb = input_data.Nb * inv2

        temp_Na[:, j1] = sum(input_data.Na)
        temp_Nb[:, j1] = sum(input_data.Nb)

        # sum_of_Na = sum_of_Na + input_data.Na

        input_data.P1 = ((sum(input_data.Na) * 1.20E-22 * input_data.Tg) / (input_data.plenum_volume))  # 6.022E23 *
        temp_P[input_data.isotope_index, j1] = input_data.P1
        # print 'first loop sum',((sum(sum_of_Na)))
        print ('for first loop', input_data.P1)


    elif j1 >= input_data.pin_failure_time and j1 < input_data.reactor_shut_down_time:

        input_data.P1 = ((sum(input_data.Na) * 1.2E-22 * input_data.Tg) / (input_data.plenum_volume))  # 6.022E23 *
        # print ('for loop two', input_data.P1)

        if input_data.P1 < input_data.P_ref:
            input_data.P1 = input_data.P_ref

        input_data.L = 3600 * (input_data.constant * (input_data.P1 - input_data.P_ref))
        input_data.L = (1 / input_data.plenum_volume) * input_data.L

        inv1 = (1 + input_data.lmbda * input_data.dt_nu + input_data.L * input_data.dt_nu) ** (-1)
        inv2 = (1 + input_data.lmbda * input_data.dt_nu + input_data.sweep_to_volume * input_data.dt_nu) ** (-1)

        input_data.Na = (input_data.rbf * input_data.dt_nu * input_data.production + input_data.Na) * inv1
        input_data.Nb = (input_data.Nb + input_data.L * input_data.Na * input_data.dt_nu) * inv2

        temp_P[input_data.isotope_index, j1] = input_data.P1
        temp_L[input_data.isotope_index, j1] = input_data.L
        temp_Na[:, j1] = sum(input_data.Na)
        temp_Nb[:, j1] = sum(input_data.Nb)


    # elif j1 >= input_data.reactor_shut_down_time or input_data.P1 > 2E6:
    elif input_data.P1 > 2E6:
        input_data.P1 = ((sum(input_data.Na) * 1.2E-22 * input_data.Tg) / (input_data.plenum_volume))  # 6.022E23 *
        # print ('for loop three', input_data.P1)

        if input_data.P1 < input_data.P_ref:
            input_data.P1 = input_data.P_ref

        input_data.L = 3600 * (input_data.constant * (input_data.P1 - input_data.P_ref))
        input_data.L = (1 / input_data.plenum_volume) * input_data.L

        inv1 = (1 + input_data.lmbda * input_data.dt_nu + input_data.L * input_data.dt_nu) ** (-1)
        inv2 = (1 + input_data.lmbda * input_data.dt_nu + input_data.sweep_to_volume * input_data.dt_nu) ** (-1)

        input_data.Na = (00 + input_data.Na) * inv1
        input_data.Nb = (input_data.Nb + input_data.L * input_data.Na * input_data.dt_nu) * inv2

        temp_P[input_data.isotope_index, j1] = input_data.P1
        temp_L[input_data.isotope_index, j1] = input_data.L
        temp_Na[:, j1] = sum(input_data.Na)
        temp_Nb[:, j1] = sum(input_data.Nb)

time_t = np.array([i4 for i4 in range(input_data.no_of_iteration)])
# print "L",input_data.L
plt.plot(input_data.dt_nu * time_t, temp_Na[input_data.isotope_index, :])
plt.show()

plt.plot(input_data.dt_nu * time_t, temp_Nb[input_data.isotope_index, :] / (4 * 1E6 * 10 * 1E10))
plt.show()

plt.plot(input_data.dt_nu * time_t, temp_P[input_data.isotope_index, :])
plt.show()

plt.plot(input_data.dt_nu * time_t, temp_L[input_data.isotope_index, :])
plt.show()


plt.plot(input_data.dt_nu*time_t,input_data.Na[13,:])
plt.show()

plt.plot(input_data.dt_nu*time_t,input_data.Nb[13,:])
plt.show()

plt.plot(input_data.dt_nu * time_t,  input_data.P)
plt.show()

plt.plot(input_data.dt_nu*time_t,input_data.L)
plt.show()




'''
