# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 09:53:14 2017

@author: Parth
"""

import math

### Bubble details

no_of_channel = 438
g = 9.8  # gravitational acceleration
H = 5  # height of the pool
kb = 1.3807E-23  # boltzman constant in J/k
Na = 6.022E23  # avogrado number
lembda = 0  # decay constant
c = 0  # total rate( diffusion, sedimentation, inertial impaction)
cp = 1.2844  # J/kg
volume_of_plenum = 0.710 * math.pi * 0.25 * (6.6e-3)**2 #m

gas_density = 5.761e3 # kg/m3 density of xenon
bubble_dia = 0.016946
vol_bubble = 0.1667 * math.pi * bubble_dia**3
bubble_mass = gas_density * vol_bubble # gas density * volume of bubble
kb = 1.3807E-23  # boltzman constant in J/k
plenum_pressure = 3e6
gas_temp = 800+273
flow_rate = 0.000158
db = []


time = 10
factor = (gas_density*no_of_channel* volume_of_plenum)
num_den_in_one_bubble = plenum_pressure * vol_bubble/ kb*gas_temp

mass_released_at_time = factor * flow_rate * time

#        num_density_released_at_time = 1

no_of_bubble_at_instance = mass_released_at_time  / bubble_mass

for i in range(int(no_of_bubble_at_instance)):
    db.append([1.69])
