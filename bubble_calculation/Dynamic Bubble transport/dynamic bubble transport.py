# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 15:37:27 2017

@author: Parth
"""


import numpy as np
import math
from matplotlib import pyplot as plt

number = int(2e4)
dt = 1E-6
samplerate = 1000

file = open('output.txt','w')



#%% Initialisations

#flow_rate = [0]
gas_temp = [1000+273]
db = []


# %%Aerosol size

# aerosol_size = np.linspace(0.01E-6,5e-6,int(30))
# aerosol_size = [10E-6]
aerosol_size = [1e-8, 1e-7, 1e-6, 3e-6]
DF = np.zeros([len(aerosol_size), number])





#%% Constants

### Temperatures and Heat Transfer Coeff
Tb = 1154.1  # Boiling point of sodium in kelvine
Ts = 800 + 273  # sodium pool temperature
T_plenum = 1000 + 273
h = 10  # check (heat transfer coefficient)

### Pressures
plenum_pressure = [5e6]  # pressure of plenum 50 bar
sodium_pressure = 1e5 # sodium ambient pressure 1 bar
P0 = 1e5  # 13# Atmosperic pressure 1 bar = 1E5 pascals

### Molar Mass
molar_mass_gas = 135e-3  # kr#135E-3 #(xe-135)
molar_mass_aerosol = 239E-3  # (fuel U-235)

### Sodium Property
rhol = 759.7575  # 949-0.223*(T[0]-273.15)-1.7E-5*(T[0]-273.15)**2 #density of sodium
mul = 0.182

### Other Things
g = 9.8  # gravitational acceleration
H = 5  # height of the pool
kb = 1.3807E-23  # boltzman constant in J/k
Na = 6.022E23  # avogrado number
lembda = 0  # decay constant
c = 0  # total rate( diffusion, sedimentation, inertial impaction)
cp = 1.2844  # J/kg
volume_of_plenum = 100e-6 #m
gas_density = 5.761e3 # kg/m3 density of xenon
cpx = 0.16# 0.16e3 # specific heat of xenon j/kg.k

#%% core number density data :--->>> Read

with open('number_density_in_core.txt') as f:
    #    count = sum(1 for line in f)
    x1 = []
    for line in f:
        one_line1 = [i for i in line.split()]
        x1.append(one_line1)

number_density_in_core = [list(list(zip(*x1))[i]) for i in range(len(x1[0]))]

number_density_in_core = number_density_in_core[0]

number_density_in_core = list(map(float, number_density_in_core))





#%% Vapor pressure data :---->>> Read

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
        partial_pressure[i] =1e5* (10 ** ((A[i]) + B[i] / T_plenum + C[i] * math.log10(T_plenum) + D[i] * T_plenum * 1e-3))
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




gas_viscosity = 25e-6 # Viscosity of xenon (as majority will be xenon)

length_of_pipe1 = 0.5 #m From failed clad location to upper plenum
length_of_pipe2 = 0.835 #m sheild location

effective_dia_of_pipe = 1.495e-3



velocity_of_gas = []

distance = [0]

s = 0

flow_rate = 0

d1 = effective_dia_of_pipe

d2 = 8.156e-3 #dia of sodium flow area

# flow_area = Triangle area - 0.5 Circle Area

flow_area1 = math.sqrt(3)/4 * (6.6e-3)**2 - 0.5 * (math.pi/4 * (6.6e-3)**2) #m2



flow_area2 = math.sqrt(3)/4 * (d2)**2 - 0.5 * (math.pi/4 * (d2)**2) #m2

volume_of_flow_path1 = flow_area1 * length_of_pipe1
volume_of_flow_path2 = flow_area2 * length_of_pipe2

file.write('Gas viscosity \t{}\n'.format(gas_viscosity))
file.write('Volume of flow path1 \t{}\n'.format(volume_of_flow_path1))
file.write('Volume of flow path2 \t{}\n'.format(volume_of_flow_path2))
file.write('flow area1 \t{}\n'.format(flow_area1))
file.write('flow area2\t{}\n'.format(flow_area2))

#%% Program Starts here

file.write('TIME\t\tFLOW RATE\t\tGAS TEMP\t\tPLENUM PRE\t\tGAS VELOCITY\t\t DISTACE\n')

for i in range(number):


    ### Time

    time = i*dt



    file.write('{:0.13f}\t'.format(time))

#    flow_rate = flow_area * effective_dia_of_pipe**2 * (plenum_pressure[i] - sodium_pressure) / (48 * length_of_pipe * gas_viscosity * volume_of_flow_path)

    ### Velocity of gas

    if distance[i] < 0.4:
        velocity_of_gas.append(flow_rate * volume_of_flow_path1 / flow_area1)
    else:
        velocity_of_gas.append(flow_rate * volume_of_flow_path2 / flow_area2)
    file.write('{}\t'.format(velocity_of_gas[i]))



    ### Distance Travelled by gas

    distance.append(s)

    file.write('{}\n'.format(distance[i]))


    '''this flow rate is basically leak rate, flow rate/ volume'''

#    flow_rate = (flow_area * (plenum_pressure[i] - sodium_pressure)/ 48)  *  (798*gas_viscosity * length_of_pipe1/d1**2 + 9 * mul * velocity_of_gas[i] * time / d2**2 )**(-1)

    file.write('{}\t'.format(flow_rate))


    ### Temperature of Gas

    gas_temp.append(gas_temp[i] + dt * (-h * flow_area1 * (gas_temp[i] - Ts))/(gas_density * volume_of_plenum * cpx))

    file.write('{}\t'.format(gas_temp[i]))


    ### Plenum Pressure

    plenum_pressure.append(plenum_pressure[i] - plenum_pressure[i] * dt * ((flow_rate)))
    #- 1 / gas_temp[i] * (-h * flow_area * (gas_temp[i] - Ts))/(gas_density * volume_of_plenum * cpx)))

    file.write('{}\t'.format(plenum_pressure[i]))

    s = s + velocity_of_gas[i] * time



    if distance[i]>0.7:
        db.append((0.976 * 6 * flow_rate**(6.0/5.0) / (g**(3.0/5.0) * math.pi) )**(1.0/3.0))
        print(db)
#        break
    else:
        pass


file.close()












#%% Plots

### Plenum Pressure
plt.figure()
plt.plot([j*dt for j in range(i)],plenum_pressure[0:-2],label = 'Plenum pressure')
plt.legend()

### Gas Temperature

'''
plt.figure()
plt.plot(gas_temp,label = 'Gas Temperature')
plt.legend()

'''
### Velocity Of Gas
plt.figure()
plt.plot([j*dt for j in range(i)],velocity_of_gas[0:-1],label = 'Gas velocity')
plt.legend()

plt.figure()
plt.plot([j*dt for j in range(i)],distance[0:-1],label = 'Distance')
plt.legend()

