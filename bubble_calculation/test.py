import numpy as np
import math
from matplotlib import pyplot as plt

number = int(1e6)
dt = 1E-4
samplerate = 1000

#%% Initialisations
T = np.zeros(number)
# T_plenum = np.zeros(number)
Db = np.zeros(number)
v = np.zeros(number)  # bubble rising velocity m/s
t = np.zeros(number)
y = np.zeros(number)
rhol = np.zeros(number)
rhog = np.zeros(number)
Pb = np.zeros(number)
N = np.zeros(number)  # chack
Ma = np.zeros(number)
RF = np.zeros(number)
Vb = np.zeros(number)
diff = np.zeros(number)
sedimentation = np.zeros(number)
inertial_impaction = np.zeros(number)
N_aero = np.zeros(number)

# %%Aerosol size
# aerosol_size = np.linspace(0.01E-6,5e-6,int(30))
# aerosol_size = [10E-6]
aerosol_size = [1e-8, 1e-7, 1e-6, 3e-6]
DF = np.zeros([len(aerosol_size), number])

#%% Constants


Tb = 1154.1  # Boiling point of sodium in kelvine
Ts = 800 + 273  # sodium pool temperature
cp = 1.2844  # J/kg
plenum_pressure = 5e6  # pressure of plenum 50 bar
g = 9.8  # gravitational acceleration
P0 = 1e5  # 13# Atmosperic pressure 1 bar = 1E5 pascals
H = 5  # height of the pool
h = 5  # check (heat transfer rate)
kb = 1.3807E-23  # boltzman constant in J/k
Na = 6.022E23  # avogrado number
lembda = 0  # decay constant
molar_mass_gas = 135e-3  # kr#135E-3 #(xe-135)
molar_mass_aerosol = 239E-3  # (fuel U-235)
c = 0  # total rate( diffusion, sedimentation, inertial impaction)
partial_pressure = 0
volume_of_plenum = 100e-6 #m
volume_of_big_bubble = 0
#%% Initial Values

T[:] = 800 + 273  # check
T_plenum = 1000 + 273
rhol[:] = 759.7575  # 949-0.223*(T[0]-273.15)-1.7E-5*(T[0]-273.15)**2 #density of sodium
Pb[0] = P0 + rhol[0] * g * (H)  # pressure of bubble
v[0] = 1e-3  # check value #velocity of bubble
# N[0] = 2.1136E27 #number density inside bubble (total = gas + aerosols)

Db[0] = 2e-2  # bubble diameter assumed
# rho_aero[0] = rhop #aerosol density
Vb[0] = math.pi * 0.1667 * Db[0] ** 3  # bubble volume

#%% Print Statements:

print('pressure is:', Pb[0], 'Pa')
print('sodium density: {:0.3e}'.format(rhol[0]), 'kg/m3')
print('bubble diameter: {:0.3e}'.format(Db[0] * 1E2), 'cm')
#%% core number density data :--->>> Read


with open('number_density_in_core.txt') as f:
    #    count = sum(1 for line in f)
    x1 = []
    for line in f:
        one_line1 = [i for i in line.split()]
        x1.append(one_line1)

number_density_in_core = [list(list(zip(*x1))[i]) for i in range(len(x1[0]))]
number_density_in_core = number_density_in_core[0]
print(number_density_in_core)
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

#partial_pressure = [10 ** ((A[i]) + B[i] / T_plenum + C[i] * math.log10(T_plenum) + D[i] * T_plenum * 1e-3) for i in range(len(Name))]

partial_pressure = np.zeros(len(Name))

for i in range(len(Name)):
    if Name[i] != 'i':
        partial_pressure[i] =1e5* (10 ** ((A[i]) + B[i] / T_plenum + C[i] * math.log10(T_plenum) + D[i] * T_plenum * 1e-3))
    else:
        partial_pressure[i] = 1e5*(10**(A[i]-B[i]/(T_plenum + C[i])))

rho_aero = np.array(rho_aero)



print('partial_pressure', partial_pressure)

pressure_of_non_condensable = Pb[0]

number_density_of_isotope = np.zeros(len(Name))
number_density_of_isotope_after_cooling = np.zeros(len(Name))

condition_volatile = []
condition_num_partial = []

for i in range(len(partial_pressure)):
    if partial_pressure[i] > plenum_pressure:
        condition_volatile.append(i)
        number_density_of_isotope[i] = number_density_in_core[i]
    else:
        condition_num_partial.append(i)
        number_density_of_isotope[i] = partial_pressure[i] * volume_of_plenum/(kb * T_plenum)

constant = sum(number_density_of_isotope)*kb*T_plenum

volume_of_big_bubble = constant / Pb[0]

diameter_of_big_bubble = (volume_of_big_bubble/(0.1667*math.pi))

no_of_bubble = 1e5

Vb[0] = volume_of_big_bubble/no_of_bubble

number_density_of_isotope = (number_density_of_isotope) / no_of_bubble

Db[0] = (sum(number_density_of_isotope) * kb * T[0] / (Pb[0] * math.pi * 0.1667))**(1.0/3.0)


#%% Partial Pressure in the sodium i.e at sodium sub-cooled temperature

#partial_pressure = [10 ** ((A[i]) + B[i] / T[0] + C[i] * math.log10(T[0]) + D[i] * T[0] * 1e-3) for i in range(len(Name))]

for i in range(len(Name)):
    if Name[i] != 'i':
        partial_pressure[i] =1e5* (10 ** ((A[i]) + B[i] / T[0] + C[i] * math.log10(T[0]) + D[i] * T_plenum * 1e-3))
    else:
        partial_pressure[i] = 1e5*(10**(A[i]-B[i]/(T[0] + C[i])))


for i in range(len(partial_pressure)):
    if i not in condition_volatile:
        number_density_of_isotope_after_cooling[i] = partial_pressure[i] * volume_of_plenum/(kb * T_plenum)

number_density_of_isotope_after_cooling = np.array(number_density_of_isotope_after_cooling)

number_density_of_isotope = np.array(number_density_of_isotope)



N_aero = number_density_of_isotope - number_density_of_isotope_after_cooling
molar_mass = np.array([83.8e-3,85.4678e-3,87.62e-3,87.62e-3,91.224e-3,101.07e-3,127.6e-3,126.9045e-3,131.293e-3,132.9055e-3,137.327e-3,138.9055e-3,140.116e-3,238.0289e-3,244e-3])
print('number_dens_of_isotope_after_cooling', number_density_of_isotope_after_cooling)
print('Number Density to be aerosol {}'.format(N_aero))
print('number_dens_of_isotope', number_density_of_isotope)
print('Volume of Small Bubble: {}'.format(Vb[0]))
print('Diameter of Small Bubble: {}'.format(Db[0]))
print ('volume of big bubble',volume_of_big_bubble)
print('big bubble diameter is : {}'.format(diameter_of_big_bubble))
#Ma[0] = (N_aero[k]*molar_mass[k])/Na









'''
for i in range(len(partial_pressure)):
    if partial_pressure[i]<Pb[0]:
        condition_to_be_aerosol.append(i)
        pressure_of_non_condensable = pressure_of_non_condensable - partial_pressure[i]
for i in range(len(Name)):
    number_density_of_isotope[i] =  partial_pressure[i]*Vb[0]/(kb*T[0])

n_non_condensable = pressure_of_non_condensable*Vb[0]/(kb*T[0])

'''
