import numpy as np
from matplotlib import pyplot as plt
import math
number = int(1e6)
#%% Initialisations
T = np.zeros(number)
Db = np.zeros(number)
A = np.zeros(number)
v = np.zeros(number) #bubble rising velocity m/s
t = np.zeros(number)
y = np.zeros(number)
rhol = np.zeros(number)
rhog = np.zeros(number)
Pb = np.zeros(number)
N = np.zeros(number) #chack
Ma = np.zeros(number)
N_aerosol = np.zeros(number)
Mg = np.zeros(number)
par = np.zeros(number)
RF = np.zeros(number)
rho_aero = np.zeros(number)
Vb = np.zeros(number)
diff = np.zeros(number)
sedimentation = np.zeros(number)
inertial_impaction = np.zeros(number)
number_dens_removed = np.zeros(number)

#%%Aerosol size
#aerosol_size = np.linspace(0.01E-6,5e-6,int(30))
#aerosol_size = [1E-6,1E-5]
aerosol_size  =[1e-8,1e-7,1e-6,3e-6]
DF = np.zeros([len(aerosol_size),number])

#%% Constants

j = 0 #second loop constant
Tb = 1154.1 # Boiling point of sodium in kelvine
Ts = 661+273 # sodium pool temperature
cp = 1284.4 #kJ/kg
k = 1.38E-23 #Boltzman constant J/K
g = 9.8 # gravitational acceleration
P0 = 1e5# Atmosperic pressure 1 bar = 1E5 pascals
H = 0#5 # height of the pool
h = 0.1 #check (heat transfer rate)
k = 1.3807E-23 #boltzman constant in J/k
Na = 6.022E23 #avogrado number
lembda = 1e-2 #decay constant
molar_mass_gas = 135E-3 #(xe-135)
molar_mass_aerosol = 239E-3 #(fuel U-235)

c = 0 # total rate( diffusion, sedimentation, inertial impaction)
mug = 2.28E-5 #gas viscosity (xenon)
#rhop =1e4 #input('enter the value of rhop:  ')


#%% Initial Values

T[0]= 27 +273 #check
rhol[0]=949-0.223*(T[0]-273.15)-1.7E-5*(T[0]-273.15)**2 #density of sodium
Pb[0] = P0+ rhol[0]*g*(H) #pressure of bubble
v[0] =1e-8 #check value #velocity of bubble
N[0] = 2.1136E27 #number density inside bubble (total = gas + aerosols)
#rhog[0] = (N[0]*molar_mass_gas/Na)
Db[0]= (1e-6/(0.1667*math.pi))**(1.0/3)#10e-2 ##1.2e-2 #bubble diameter
#rho_aero[0] = rhop #aerosol density
Vb[0] = math.pi * 0.1667 * Db[0]**3 #bubble volume
Ma[0] = rho_aero[0]*Vb[0] #initial mass of aerosols in bubble
#%%

with open('vapor_pressure.txt') as f:
#    count = sum(1 for l in f)
    x = []
    for line in f:
        y = [i for i in line.split()]
        x.append( y )
f.close()

Name,A,B,C,D = [list(list(zip(*x))[i]) for i in range(5)]

A = list(map(float,A))
B = list(map(float,B))
C = list(map(float,C))
D = list(map(float,D))
#T = [800 for i in range(len(Name))]
#partial_pressure = []

partial_pressure=[10**((A[i]) + B[i]/T[0] + C[i]*math.log10(T[0]) + D[i]* T[0]*1e-3) for i in range(len(Name))]

partial_pressure = np.asarray(partial_pressure)

pressure_of_non_condensable = Pb[0] - sum(partial_pressure[:])

partial_pressure = 1e5*partial_pressure


molar_mass = [85,132,87,137,101,239,235,80,84,131,107,56]


n0,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11 = [partial_pressure[i]*Vb[0]/(k*T[0]) for i in range(len(Name))]

n_non_condensable = pressure_of_non_condensable*Vb[0]/(k*T[0])

N[0] = n0+n1+n2+n3+n4+n5+n6+n7+n8+n9+n10+n11+n_non_condensable

rho_aero[0] = (n5+n6)*molar_mass_aerosol/Na

num_de = [n0,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11]




j = 0
print('Bubble Diameter:{}'.format(Db[0]),'m \n')
print('Pressure:{:0.3e}'.format(Pb[0]),'Pa \n')
print('Bubble Volume:{}'.format(Vb[0]),'m3 \n')




for i in num_de:
    print('{}:{:0.1e}'.format(Name[j],i))
    j = j+1




j = 0





file = open('number_density.txt','w')
file.write('#Bubble Diameter:{}'.format(Db[0])+' m \n')
file.write('#Pressure:{:0.3e}'.format(Pb[0])+' Pa \n')
file.write('#Bubble Volume:{}'.format(Vb[0])+' m3 \n')
for i in num_de:
    file.write('{}:{:0.1e}\n'.format(Name[j],i))
    j = j+1
file.close()
