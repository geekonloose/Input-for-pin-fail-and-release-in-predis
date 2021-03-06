# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 15:37:27 2017

@author: Parth
"""
# TODO: Implement the destruction condition and condition to reach on top
# formula used to calculate partial pressure
    # (10 ** ((A[i])+ B[i] / T_plenum + C[i] * math.log10(T_plenum) + D[i] * T_plenum * 1e-3)



import numpy as np
import math
from matplotlib import pyplot as plt
plt.close('all')

error = 0
number = int(1e2)
dt = 1e-2
samplerate = 1000

# file = open('output.txt','w')


# %% Initialisations
volume_of_plenum = np.zeros(number)
# flow_rate = [0]
gas_temp = [1000+273]
db = []
test_flow = []
per = []
number_density_released = []
dia = []
no_of_bubble_generated = []


Ma = []
c = []
RF = []
diff = []
sedimentation = []
inertial_impaction = []
dp = 1e-5
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
h = 1e2  # check (heat transfer coefficient)

###  Pressures
plenum_pressure = [3e6]  # pressure of plenum 50 bar
sodium_pressure = 1e5 # sodium ambient pressure 1 bar
P0 = 1e5  # 13# Atmosperic pressure 1 bar = 1E5 pascals

### Molar Mass
molar_mass_gas = 135e-3  # kr#135E-3 #(xe-135)
molar_mass_aerosol = 239E-3  # (fuel U-235)

### Sodium Property
rhol = 759.7575 #sodium density at 800 C
# 949-0.223*(T[0]-273.15)-1.7E-5*(T[0]-273.15)**2 #density of sodium
mul = 0.182e3

### Other Things
no_of_channel = 438
g = 9.8   # gravitational acceleration
H = 5   # height of the pool
kb = 1.3807E-23   # boltzman constant in J/k
Na = 6.022E23   # avogrado number
lembda = 0   # decay constant
#c = 0   # total rate( diffusion, sedimentation, inertial impaction)
cp = 1.2844   # J/kg
volume_of_plenum[0] = 0.710 * math.pi * (6.6e-3)**2  # m
gas_density = 5.761e3   # kg/m3 density of xenon
cpx = 0.16   # 0.16e3  # specific heat of xenon j/kg.k
mug = 25e-6
### Bubble details


bubble_dia = 0.016946
vol_bubble = 0.1667 * math.pi * bubble_dia**3
bubble_mass = gas_density * vol_bubble   # gas density * volume of bubble






# %% NU
#flow_rate = 1e-4
#
#
#nu = flow_rate / volume_of_plenum[0]   # release from pin\

nu = 1e-1

#print(nu*3e6)

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
        partial_pressure[i] = (1e5 * (10 ** ((A[i])+ B[i] / T_plenum + C[i] * math.log10(T_plenum) + D[i] * T_plenum * 1e-3)))
    else:
        partial_pressure[i] = 1e5*(10**(A[i]-B[i]/(T_plenum + C[i])))


### Numbeer Density in plenum

condition_volatile = []
condition_num_partial = []

number_density_of_isotope = [[0 for i in range(len(Name))] for j in range(1)]


for i in range(len(partial_pressure)):

    if partial_pressure[i] > plenum_pressure:

        condition_volatile.append(i)

        number_density_of_isotope[0][i] = number_density_in_core[i]

    else:

        condition_num_partial.append(i)

        number_density_of_isotope[0][i] = float(partial_pressure[i] * volume_of_plenum[0]/(kb * T_plenum))

#constant = sum(number_density_of_isotope)*kb*T_plenum

plenum_number_density = sum(number_density_of_isotope[0])



### Pre Data For Flow Calculation




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

# file.write('Gas viscosity \t{}\n'.format(gas_viscosity))
# file.write('Volume of flow path1 \t{}\n'.format(volume_of_flow_path1))
# file.write('Volume of flow path2 \t{}\n'.format(volume_of_flow_path2))
# file.write('flow area1 \t{}\n'.format(flow_area1))
# file.write('flow area2\t{}\n'.format(flow_area2))






num_den_remained = []


# %%
# %%
#file = open('num_den_released.txt','w')
# %% Program Starts here

 #file.write('TIME\t\tFLOW RATE\tvelocity\tDISTACE\tTEMP\tPRESSURE\n')
k = 0
j = 0


molar_mass = np.array([83.8e-3,85.4678e-3,87.62e-3,87.62e-3,91.224e-3,101.07e-3,127.6e-3,126.9045e-3,131.293e-3,132.9055e-3,137.327e-3,138.9055e-3,140.116e-3,238.0289e-3,244e-3])




initial_mass = []

for i in range(len(Name)):

    initial_mass.append(number_density_of_isotope[0][i] * molar_mass[i] / Na)

#====================================================================================================
#====================================================================================================
###                                    PPROGRAM STARTS HERE
#====================================================================================================
#====================================================================================================

for i in range(number):



### Time

    time = i*dt

### Flow rate

    "HERE RELEASE RATE = FLOW RATE (M3/S) / VOLUME OF PLENUM"

#    number_density_of_isotope[i] = number_density_of_isotope[i] - number_density_released[i]

    number_density_released.append(list(np.array(number_density_of_isotope[i])*dt*nu))


    num_den_remained.append(list(np.array( number_density_of_isotope[i])- np.array(number_density_released[i])))



    dia.append((6 * sum(number_density_released[i]) * kb * T_plenum / (math.pi * plenum_pressure[i]))**(1/3))


    if dia[i] > 8e-3:
        no_of_bubble_generated.append(int(dia[i] / d1))
    else:
        no_of_bubble_generated.append(1)


    for i2 in range(int(no_of_bubble_generated[i])):
            db.append([0.016946])
            Ma.append([[i] for i in initial_mass])
            RF.append([[] for i in initial_mass])
            diff.append([[0 for i11 in range(i)] for i12 in range(len(Name))])
            sedimentation.append([[0 for i11 in range(i)] for i12 in range(len(Name))])    #change here
            inertial_impaction.append([[0 for i11 in range(i)] for i12 in range(len(Name))])
            c.append([[0 for i11 in range(i)] for i12 in range(len(Name))])



### Gas temp, plenum pressure, vol of plenum

    plenum_pressure.append(3e6*math.exp(-nu*time))

    volume_of_plenum[i] = sum(num_den_remained[i]) * kb * T_plenum / (plenum_pressure[i])

    gas_temp.append(gas_temp[0] * math.exp(-2*nu*time))




### Partial pressure calculation at each steps

    for i1 in range(len(Name)):

        if Name[i1] != 'i':
            partial_pressure[i1] = (1e5 * (10 ** ((A[i1])+ B[i1] / gas_temp[i] + C[i1] * math.log10(gas_temp[i]) + D[i1] * gas_temp[i] * 1e-3)))
        else:
            partial_pressure[i1] = 1e5*(10**(A[i1]-B[i1]/(gas_temp[i] + C[i1])))

### Number Density calculation at each steps

    condition_volatile = []
    condition_num_partial = []

    n1 = []

    for j1 in range(len(partial_pressure)):

        if partial_pressure[j1] > plenum_pressure[i]:

            condition_volatile.append(j1)
            n1.append( number_density_in_core[j1] )

        else:

            condition_num_partial.append(j1)

            n1.append(np.array(partial_pressure[j1] * volume_of_plenum[i]/(kb * T_plenum)) - (np.array(number_density_released[i][j1])))


    number_density_of_isotope.append(n1)


### bubble transport starts here


#%% Bubble Module

### Fix for equalling the size of db for late bubbles

    for i6 in range(1,no_of_bubble_generated[i]+1):

        for i5 in range(len(db[0])-1):

            db[-i6].insert(0,0)


    if i == 0:    # to make sure that variable defined only once

        y = [[0] for j2 in range(len(db))]

        v = [[0.001 for i4 in range(1)] for j2 in range(len(db))]

        Pb = [[1e5 for i4 in range(1)] for j2 in range(len(db))]

        rhog = [[0 for i4 in range(1)] for j2 in range(len(db))]


### Fix for difference in length of list

    while len(db) != len(y):
        y.append([0 for i in range(len(y[0]))])
        v.append([0.001 for i in range(len(y[0]))])
        Pb.append([1e5 for i in range(len(y[0]))])
        rhog.append([0 for i in range(len(y[0]))])





#====================================================================================================
#====================================================================================================
####                                     FOR ALL BUBBLES
#===================================================================================================
#====================================================================================================



    for i3 in range(len(db)):    #for each bubble diameter

        Pb[i3].append( P0 + rhol * g * ( H - y[i3][i]) )

        if i > 0 :

            db[i3].append( ( 6 * sum( number_density_of_isotope[i] ) * kb * time / ( 3.14 * Pb[i3][i] ) )**( 1 / 3.0 )) # diameter equation
#            print('i3',i3,'i',i,'db',db)

        rhog[i3].append(sum(number_density_of_isotope[i]) * molar_mass_gas / Na)

        v[i3].append( v[i3][i] + dt * ( 6 / ( rhol * 3.14 * db[i3][i]**3 ) ) * ( ( 3.14 * db[i3][i]**3 * g * ( rhol - rhog[i3][i] ) ) / 6 - 12 * 3.14 * db[i3][i] * v[i3][i] ) ) #velocity of bubble

        y[i3].append(i * dt * v[i3][i])


#====================================================================================================
#====================================================================================================
###                               AEROSOL DEPOSITION MODULE
#====================================================================================================
#====================================================================================================

### Initialisation of the aerosol mass


        for i9 in range(len(Name)):




            lembd = kb * gas_temp[i] / ( float(math.sqrt(2) ) * math.pi * Pb[i3][i] * db[i3][i]**2 )    #mean free path of a gas molecule

            cunn = 1 + ( 2 * lembd / dp ) * ( 1.257 + 0.4 * math.exp( - 0.55 * dp / lembd ) ) #cunningham slip correction

            theta = kb * gas_temp[i] * cunn / ( 3 * 3.14 * mug * dp )

            tau = rho_aero[i9] * dp**2 * cunn / ( 18 * mug )

            diff[i3][i9].append( 1.8 * ( 8 * theta / ( v[i3][i] * db[i3][i]**3 ) )**( 1 / 2.0 ) )    #diffusion coefficient

            sedimentation[i3][i9].append( 1.5 * g * tau / ( db[i3][i] * v[i3][i] ) )    #sedimentation coefficient

            inertial_impaction[i3][i9].append( 18 * v[i3][i] * tau / db[i3][i]**2 )    # inertial impaction coefficient


            c[i3][i9].append( diff[i3][i9][i] + sedimentation[i3][i9][i] + inertial_impaction[i3][i9][i] )


            Ma[i3][i9].append( Ma[i3][i9][0] * math.exp( - y[i3][i] * c[i3][i9][i] ))    #calculates aerosol masss

            RF[i3][i9].append( math.exp( - y[i3][i] * c[i3][i9][i] ) )


#%% Plots



plt.figure()
plt.plot(np.linspace(1., number, number ) * dt,RF[10][4][:],label = 'release fraction')
plt.xlabel('Time')
plt.ylabel('Release fracton (Mass in bubble/ Total mass)')



plt.figure()
plt.plot(np.linspace(1., number, number ) * dt,np.array(db[10][:]) * 100,label = 'release fraction')
plt.xlabel('Time')
plt.ylabel('dia')








'''


#if error == 0:
#a = int(input('enter'))
#if a==1:
plt.plot([dt*j for j in range(i+2)],plenum_pressure,label = "plenum pressure")
plt.legend()

### plot Gas Temperature
plt.figure()
plt.plot(gas_temp,label = 'Gas Temperature')
plt.legend()


plt.figure()
plt.plot(volume_of_plenum,label = "plenum vol")
plt.legend()

'''





'''
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
'''
