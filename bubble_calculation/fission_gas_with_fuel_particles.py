#-*- coding: utf-8 -*-
"""
Created on Sat Jul 29 17:10:21 2017
@author: parth
"""

import numpy as np
import math
from matplotlib import pyplot as plt

number = int(1e6)
dt = 1E-4
samplerate = 1000

#%% Initialisations
T = np.zeros(number)
Db = np.zeros(number)
v = np.zeros(number) #bubble rising velocity m/s
t = np.zeros(number)
y = np.zeros(number)
rhol = np.zeros(number)
rhog = np.zeros(number)
Pb = np.zeros(number)
N = np.zeros(number) #chack
Ma = np.zeros(number)
Mg = np.zeros(number)
par = np.zeros(number)
RF = np.zeros(number)
rho_aero = np.zeros(number)
Vb = np.zeros(number)
diff = np.zeros(number)
sedimentation = np.zeros(number)
inertial_impaction = np.zeros(number)
number_dens_removed = np.zeros(number)
N_aero = np.zeros(number)
#%%Aerosol size
#aerosol_size = np.linspace(0.01E-6,5e-6,int(30))
#aerosol_size = [10E-6]
aerosol_size  =[1e-8,1e-7,1e-6,3e-6]
DF = np.zeros([len(aerosol_size),number])

#%% Constants


Tb = 1154.1 # Boiling point of sodium in kelvine
Ts = 800+273 # sodium pool temperature
cp = 1.2844 #J/kg
k = 1.38E-23 #Boltzman constant J/K
g = 9.8 # gravitational acceleration
P0 =1e5#13# Atmosperic pressure 1 bar = 1E5 pascals
H = 5 # height of the pool
h = 5 #check (heat transfer rate)
k = 1.3807E-23 #boltzman constant in J/k
Na = 6.022E23 #avogrado number
lembda = 1e-10 #decay constant
molar_mass_gas = 135e-3# kr#135E-3 #(xe-135)
molar_mass_aerosol = 239E-3 #(fuel U-235)
c = 0 # total rate( diffusion, sedimentation, inertial impaction)
#mug =
#mug = 2.28E-5 #gas viscosity (xenon)
#rhop =1e4 #input('enter the value of rhop:  ')
partial_pressure = 0

#%% Initial Values

T[:]= 800+273 #check
rhol[:]= 759.7575#949-0.223*(T[0]-273.15)-1.7E-5*(T[0]-273.15)**2 #density of sodium
Pb[0] = P0+ rhol[0]*g*(H) #pressure of bubble
v[0] =1e-3 #check value #velocity of bubble
#N[0] = 2.1136E27 #number density inside bubble (total = gas + aerosols)

Db[0]= 5e-2 #bubble diameter assumed
#rho_aero[0] = rhop #aerosol density
Vb[0] = math.pi * 0.1667 * Db[0]**3 #bubble volume


#%% Print Statements:

print('pressure is:',Pb[0],'Pa')
print('sodium density: {:0.3e}'.format(rhol[0]),'kg/m3')
print('bubble diameter: {:0.3e}'.format(Db[0]*1E2),'cm')





#%% Vapor pressure data

with open('vapor_pressure.txt') as f:
#    count = sum(1 for line in f)
    x = []
    for line in f:
        one_line = [i for i in line.split()]
        x.append(one_line )

Name,A,B,C,D = [list(list(zip(*x))[i]) for i in range(5)]
A = list(map(float,A))
B = list(map(float,B))
C = list(map(float,C))
D = list(map(float,D))
partial_pressure=[10**((A[i]) + B[i] / T[0] + C[i]*math.log10(T[0]) + D[i]* T[0]*1e-3) for i in range(len(Name))]

partial_pressure = np.asarray(partial_pressure)
partial_pressure = 1e5*partial_pressure
pressure_of_non_condensable = Pb[0]

for i in range(len(partial_pressure)):

    if partial_pressure[i]<Pb[0]:
        print(i)
        pressure_of_non_condensable = pressure_of_non_condensable - partial_pressure[i]






n0,n1,n2,n3,n4,n5,n6,n7 = [partial_pressure[i]*Vb[0]/(k*T[0]) for i in range(len(Name))]

n_non_condensable = pressure_of_non_condensable*Vb[0]/(k*T[0])



#%% Mass of aerosol

molar_mass = [132e-3,85e-3,86e-3,88e-3,84e-3,131e-3,235e-3,239e-3]

N[:] = n_non_condensable
N_aero[0] = n6+n7+n2+n3
rho_aero[:] = 1.2e4


Ma[0] = (n2*molar_mass[2] + n3*molar_mass[3] + n6*molar_mass[6] + n7*molar_mass[7])/Na #1e2#(N_aero[0])*molar_mass_aerosol/Na

#+ n3*molar_mass[3]


#239e-3*N_aero[0]/Na







#Ma[0] = 1e-2 #initial mass of aerosols in bubble
rhog[0] = (N[0]*molar_mass_gas/Na)
print('Aerosol Density:{:0.1e}'.format(rho_aero[0]),'kg/m3')
#mu_gas = [0.67e-3, 0.68e-3, 1.59e-3,1.86e-3,6.1e-3,6e-3,6.5e-3,0.00095,1.2e-2,0,0]
'''
s1 = 0
s2 = 0
for i in range(len(Name)):
    s1 = s1 + partial_pressure[i]*mu_gas[i]*molar_mass[i]**0.5
    s2 = s2 + partial_pressure[i]*molar_mass[i]**0.5
s1 = s1 + pressure_of_non_condensable*mu_gas[i]*molar_mass[i]**0.5
s2 = s2 + pressure_of_non_condensable*molar_mass[i]**0.5
mug = s1/s2
'''
mug = 25e-6

#%% file write

j = 0

num_de = [n0,n1,n2,n3,n4,n5,n6,n7]
file = open('number_density.txt','w')
file.write('#Bubble Diameter:{}'.format(Db[0])+' m \n')
file.write('#Pressure:{:0.3e}'.format(Pb[0])+' Pa \n')
file.write('#Bubble Volume:{}'.format(Vb[0])+' m3 \n')
for i in num_de:
    file.write('{}:{:0.1e}\n'.format(Name[j],i))
    j = j+1
file.close()

#%%  Program loop
i = 0 #initialize
j = 0 #initialize
for dp in aerosol_size:
    print('aerosol size is:',dp,'m')
#    rho_aero[:] = Ma[0]/(0.1667*math.pi*dp**3)
#    print('aerosol_density',rho_aero[0])
    for i in range(number-3):

        t[i] = i*dt

        y[i] = (i)*dt*v[i]

        N[i+1] = N[i] - lembda*dt*N[i]

        if i>0:
            Ma[i] = Ma[0]*math.exp(-y[i]*c)

        rhol[i] = 949-0.223*(T[i]-273.15)-1.7E-5*(T[i]-273.15)**2 #density of sodium

        Pb[i] =P0 + rhol[i]*g*(H-y[i]) # pressure equation

        if i > 0 :
            Db[i] = (6*N[i]*k*T[i]/(3.14*Pb[i]))**(1/3.0) # diameter equation

        Vb[i] = math.pi*0.1667*Db[i]**3

#        rho_aero[i] = Ma[i]/(0.1667*math.pi*dp**3)

        rhog[i] =  (N[i]*molar_mass_gas/Na)

#        number_dens_removed[i] = (- (rho_aero[i]-rho_aero[i-1])*Na /   molar_mass_aerosol)

        if i>0:
            number_dens_removed[i] = (Ma[i]-Ma[i-1])*Na / 0.235

        N_aero[i+1] = N_aero[i] + number_dens_removed[i] *  dt

#        T[i+1] = T[i] + dt* (0-h*3.14*Db[i]**2 * (T[i]-Ts))/(rhog[i]*Vb[i]*cp) # temperature equation


        v[i+1] = v[i]+ dt*(6/(rhol[i]*3.14*Db[i]**3))*(  (3.14*Db[i]**3 *g*(rhol[i]-rhog[i]))/6 - 12*3.14*Db[i]*v[i]) #velocity of bubble

        lembd = k*T[i]/(float(math.sqrt(2))*math.pi*Pb[i]*Db[i]**2) #mean free path of a gas molecule

        cunn = 1 + (2*lembd/dp) * (1.257 + 0.4 * math.exp(-0.55*dp/lembd)) #cunningham slip correction

        theta = k*T[i]*cunn/(3*3.14*mug*dp)

        tau = rho_aero[i]*dp**2 *cunn /(18*mug)

        diff[i] = 1.8*(8*theta/(v[i]*Db[i]**3))**(1/2.0) #diffusion coefficient

        sedimentation[i] = 1.5*g*tau/(Db[i]*v[i]) #sedimentation coefficient

        inertial_impaction[i] = 18*v[i]*tau/Db[i]**2 # inertial impaction coefficient

        c = diff[i] + sedimentation[i] + inertial_impaction[i]

        dy = dt*v[i]
        par[i] = dy*(-c*Ma[i])
        Ma[i+1] = Ma[i] + par[i]

        RF[i] = math.exp(-y[i]*c)

        if y[i]>H:
            print('bubble reached the pool surface')
            break
        elif N[i]<0:
            print('bubble gas decayed out')
            break

    DF[j,0:i] = RF[0:i]

    j = j+1


#%% Plots starts here

for j in range(len(aerosol_size)):
    plt.plot(t[0:i-1:samplerate],DF[j,0:i-1:samplerate],label='{:0.1e}'.format(aerosol_size[j]))
plt.legend()
plt.xlabel('time(second)',color='blue',size=14)
plt.ylabel('fraction in bubble',color='blue',size=14)
plt.savefig('fraction of mass in bubble')
#plt.close()
'''
plt.figure()
plt.plot(aerosol_size,1/DF[:,i-1])
plt.semilogx()
plt.xlabel('aerosol size(micro meter)')
plt.ylabel('DF')
plt.savefig('DF factor')
'''
'''
plt.figure()
plt.plot(aerosol_size,DF[:,i-1])
plt.semilogx()
plt.semilogy()

plt.plot(aerosol_size*1e6,DF[:j])
plt.semilogx()
plt.semilogy()
'''
plt.figure()
plt.plot(t[0:i],Db[0:i],label='diameter',color='blue')
plt.xlabel('time(seconds)',color='blue',size=14)
plt.ylabel('diameter (meter)',color='blue',size=14)
plt.legend()
plt.savefig('diameter')

#plt.figure()
'''
plt.plot([i for i in range(0,j)],(DF[0:j]),label='release fraction')

plt.xlabel('Time',color='blue',size=14)
plt.ylabel('fraction in bubble',color='blue',size=14)
plt.legend()
plt.savefig('Release')
'''

plt.figure()
plt.plot(t[:i],T[:i],label = 'Temperature',color='blue')
plt.xlabel('time',color='blue',size=14)
plt.ylabel('Temp (K)',color='blue',size=14)
plt.legend()
plt.savefig('temperature')
plt.figure()
plt.plot(t[:i],Pb[:i],label='pressure',color='blue')
plt.xlabel('time(seconds)',color='blue',size=14)
plt.ylabel('Pressure (Pa)',color='blue',size=14)
plt.legend()
plt.savefig('pressure')
plt.figure()
plt.plot(t[:i],v[:i],label='velocity',color='blue')
plt.xlabel('time(second)',color='blue',size=14)
plt.ylabel('velocity(m/s)',color='blue',size=14)

plt.legend()
plt.savefig('velocity')

'''
plt.figure()
plt.plot(t[:i],N[:i],label='number density',color='blue')
plt.xlabel('time(second)',color='blue',size=14)
plt.ylabel('number density',color='blue',size=14)

plt.legend()
plt.savefig('number density')
'''

#convert meter to cm
Db =(Db*1.2E2)
y = y*1.2E2-250
#convert array to list
Db = list(Db)
y = list(y)
import turtle

def draw_circle(turtle, color, size, x, y):
    turtle.penup()
    turtle.color(color)

    turtle.goto(x,y)
    turtle.pendown()

    turtle.circle(size)


turtle.shape("circle")
#turtle.resizemode("user")
turtle.shapesize(1E-2,1E-2)
#tommy.speed(7)
scale = 5000
# Draw three circles:
for i in range(int(len(y)/scale)):
    draw_circle(turtle, "green", Db[scale*i], 25, y[scale*i])





