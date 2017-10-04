# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 15:37:27 2017

@author: Parth
"""


import numpy as np
import math
from matplotlib import pyplot as plt
plt.close('all')


number = int(1e6)
dt = 1e-3
samplerate = 1000

# file = open('output.txt','w')


# %% Initialisations

# flow_rate = [0]
gas_temp = [1000+273]
db = []
test_flow = []
per = []

# %% Aerosol size

# aerosol_size = np.linspace(0.01E-6,5e-6,int(30))
# aerosol_size = [10E-6]
aerosol_size = [1e-8, 1e-7, 1e-6, 3e-6]
# DF = np.zeros([len(aerosol_size), number])


# %% Constants

### Temperatures and Heat Transfer Coeff
Tb = 1154.1  # Boiling point of sodium in kelvine
Ts = 800 + 273  # sodium pool temperature
T_plenum = 1000 + 273
h = 1e3  # check (heat transfer coefficient)

###  Pressures
plenum_pressure = [3e6]  # pressure of plenum 50 bar
sodium_pressure = 1e5 # sodium ambient pressure 1 bar
P0 = 1e5  # 13# Atmosperic pressure 1 bar = 1E5 pascals

### Molar Mass
molar_mass_gas = 135e-3  # kr#135E-3 #(xe-135)
molar_mass_aerosol = 239E-3  # (fuel U-235)

### Sodium Property
rhol = 759.7575
# 949-0.223*(T[0]-273.15)-1.7E-5*(T[0]-273.15)**2 #density of sodium
mul = 0.182e3

### Other Things
no_of_channel = 438
g = 9.8   # gravitational acceleration
H = 5   # height of the pool
kb = 1.3807E-23   # boltzman constant in J/k
Na = 6.022E23   # avogrado number
lembda = 0   # decay constant
c = 0   # total rate( diffusion, sedimentation, inertial impaction)
cp = 1.2844   # J/kg
volume_of_plenum = 0.710 * math.pi * 0.25 * (6.6e-3)**2  # m
gas_density = 5.761e3   # kg/m3 density of xenon
cpx = 0.16   # 0.16e3  # specific heat of xenon j/kg.k

### Bubble details


bubble_dia = 0.016946
vol_bubble = 0.1667 * math.pi * bubble_dia**3
bubble_mass = gas_density * vol_bubble   # gas density * volume of bubble


# %% core number density data :--->>> Read

with open('number_density_in_core.txt') as f:
    #    count = sum(1 for line in f)
    x1 = []
    for line in f:
        one_line1 = [i for i in line.split()]
        x1.append(one_line1)

number_density_in_core = [list(list(zip(*x1))[i]) for i in range(len(x1[0]))]

number_density_in_core = number_density_in_core[0]

number_density_in_core = list(map(float, number_density_in_core))





# %% Vapor pressure data :---->>> Read

with open('vapor_pressure.txt') as f:
    #    count = sum(1 for line in f)
    x = []
    for line in f:
        one_line = [i for i in line.split()]
        x.append(one_line)

Name, rho_aero, A, B, C, D = [list(list(zip(*x))[i]) for i in range(len(x[0]))]

A = list(map(float, A))
B = list(map(float, B))
C = list(map(float, C))
D = list(map(float, D))
rho_aero = list(map(float, rho_aero))

partial_pressure = np.zeros(len(Name))

for i in range(len(Name)):

    if Name[i] != 'i':
        partial_pressure[i] = (1e5* (10 ** ((A[i])+ B[i] / T_plenum + C[i] * math.log10(T_plenum) + D[i] * T_plenum * 1e-3)))
    else:
        partial_pressure[i] = 1e5*(10**(A[i]-B[i]/(T_plenum + C[i])))



#%% Numbeer Density in plenum

condition_volatile = []
condition_num_partial = []

number_density_of_isotope = [[0 for i in range(len(Name))] for j in range(1)]


for i in range(len(partial_pressure)):

    if partial_pressure[i] > plenum_pressure:
        condition_volatile.append(i)

        number_density_of_isotope[0][i] = number_density_in_core[i]

    else:

        condition_num_partial.append(i)

        number_density_of_isotope[0][i] = partial_pressure[i] * volume_of_plenum/(kb * T_plenum)

#constant = sum(number_density_of_isotope)*kb*T_plenum

plenum_number_density = sum(number_density_of_isotope[0])





#%% Pre Data For Flow Calculation




gas_viscosity = 5e-5 # Viscosity of xenon (as majority will be xenon)

l1 = 0.530 #m From failed clad location to upper plenum
l2 = 0.835 #m sheild location

effective_dia_of_pipe = 1.495e-3



velocity_of_gas = [0]

distance = [0]

s = 0

flow_rate = 0

d1 = effective_dia_of_pipe

d2 = 8.156e-3 #dia of sodium flow area

# flow_area = Triangle area - 0.5 Circle Area

# Hexagon area = 113.997e-4 m2

### Area of Flow

flow_area1 = 113.997e-4 - 217 * math.pi /4 * (6.6e-3)**2

flow_area2 = 113.997e-4 - 7 * math.pi / 4 * (36e-3)**2

### Volume of flow area

volume_of_flow_path1 = flow_area1 * l1
volume_of_flow_path2 = flow_area2 * l2

### equivalent inertia length

total_ineria_len = l1/flow_area1 + l2/flow_area2

 #file.write('Gas viscosity \t{}\n'.format(gas_viscosity))
 #file.write('Volume of flow path1 \t{}\n'.format(volume_of_flow_path1))
 #file.write('Volume of flow path2 \t{}\n'.format(volume_of_flow_path2))
 #file.write('flow area1 \t{}\n'.format(flow_area1))
 #file.write('flow area2\t{}\n'.format(flow_area2))



#%% Percentage of fission gas in the plenum

for i in range(len(number_density_of_isotope[0])):
    per.append(number_density_of_isotope[0][i]*100 / sum(number_density_of_isotope[0] ))









#%% Program Starts here

 #file.write('TIME\t\tFLOW RATE\tvelocity\tDISTACE\tTEMP\tPRESSURE\n')
k = 0
j = 0
for i in range(number):

### Time

    time = i*dt



     #file.write('{:0.10f}\t'.format(time))



### Flow rate

    '''this flow rate is basically leak rate, flow rate/ volume'''


    if distance[i]<0.8:
        sw = 0.8-distance[i]
    else:
        sw = 0


    if time<= 0 :
        time = 1
    flow_rate = (flow_rate + dt * 1 / total_ineria_len * (-gas_density * g * (l1+l2) + (flow_rate**2 / 2*gas_density) * (1 / flow_area1**2 - 1 / flow_area2**2 ) + (plenum_pressure[i]-sodium_pressure ) - 12 * gas_viscosity*gas_viscosity*l1/d1 -(12 * gas_viscosity * l1*flow_rate * distance[i]/(time * d1))))

    flow_rate = flow_rate * (1 / (gas_density*no_of_channel* volume_of_plenum))





#    flow_rate = ((math.pi/4 * d2**2 )* (plenum_pressure[i] - sodium_pressure)/ 48)  *  (no_of_channel*gas_viscosity * l1/d1**2  +9* mul * sw / d2**2 )**(-1)

#    flow_rate = flow_rate / (no_of_channel*volume_of_plenum)

#    flow_rate = flow_area1 * effective_dia_of_pipe**2 * (plenum_pressure[i] - sodium_pressure) / (48 * l1 * gas_viscosity * volume_of_flow_path1)

#    flow_rate = flow_rate / (no_of_channel*volume_of_plenum)

    test_flow.append(flow_rate)

#    print('{:0.10e}'.format(flow_rate))

     #file.write('{:0.2e}\t'.format(flow_rate))





### Velocity of gas

    if distance[i] < 0.4:
        velocity_of_gas.append(flow_rate * volume_of_flow_path1 / flow_area1)
    else:
        velocity_of_gas.append(flow_rate * volume_of_flow_path2 / flow_area2)

     #file.write('{:0.2f}\t'.format(velocity_of_gas[i]))





### Distance Travelled by gas

    distance.append(s)

     #file.write('{:0.2f}\t'.format(distance[i]))




### Temperature of Gas

    gas_temp.append(gas_temp[i] + dt * (-h * flow_area2 * (gas_temp[i] - Ts))/(gas_density * volume_of_plenum * cpx))

     #file.write('{:0.2f}\t'.format(gas_temp[i]))


### Plenum Pressure

    plenum_pressure.append(plenum_pressure[i] - plenum_pressure[i] * dt * ((flow_rate) -1 / gas_temp[i] * (gas_temp[i+1]-gas_temp[i])))

     #file.write('{:0.2f}\n'.format(plenum_pressure[i]))

    s = s + velocity_of_gas[i] * time








### condition to stop

    if distance[i] > 1.5 and distance[i] < 1.5 + 0.8:

        print('Bubble reached top of assembly')

        factor = (gas_density*no_of_channel* volume_of_plenum)

#        factor = 1

#        db.append((0.976 * 6 * (flow_rate * factor )**(6.0/5.0) / (g**(3.0/5.0) * math.pi) )**(1.0/3.0))



    '''

### Number of bubble generated at instance

        num_den_in_one_bubble = plenum_pressure[i] * vol_bubble/ kb*gas_temp[i]

        mass_released_at_time = factor * flow_rate * time

#        num_density_released_at_time = 1

        no_of_bubble_at_instance = mass_released_at_time  / bubble_mass

        for i in range(int(no_of_bubble_at_instance)):
            db.append([0.016946])
    '''



### change of number density as pressure changes

    n1 = []

    for j in range(len(partial_pressure)):

        if partial_pressure[j] > plenum_pressure[i]:
            condition_volatile.append(j)
            n1.append( number_density_in_core[j] )

        else:

            condition_num_partial.append(j)

            n1.append( partial_pressure[j] * volume_of_plenum/(kb * T_plenum))


    number_density_of_isotope.append(n1)


#        remaining_num_density = sum(number_density_of_isotope[-1]) -





    if distance[i] > 1.5 + 0.8:
        break
#print('Bubble dia meter:',db,'m')

 #file.close()





#import os
#os.start #file('output.txt')







#%% Plots


### plot flow rate

plt.figure()
plt.plot(test_flow,label = 'flow rate')
plt.legend()
plt.savefig('flow rate')


### plot Plenum Pressure
plt.figure()
plt.plot([j*dt for j in range(i)],plenum_pressure[0:-2],label = 'Plenum pressure')
plt.legend()
plt.savefig('plenum pressure')


### plot Gas Temperature
plt.figure()
plt.plot(gas_temp,label = 'Gas Temperature')
plt.legend()
plt.savefig('gas temperature')


### plot Velocity Of Gas
plt.figure()
plt.plot([j*dt for j in range(i)],velocity_of_gas[0:-2],label = 'Gas velocity')
plt.legend()
plt.savefig('velocity of gas')


### plot Distance
plt.figure()
plt.plot([j*dt for j in range(i)],distance[0:-2],label = 'Distance')
plt.legend()

plt.savefig('distance')

