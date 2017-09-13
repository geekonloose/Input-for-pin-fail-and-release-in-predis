
import numpy as np
import math
from matplotlib import pyplot
'''
# Input data
R = 6.6  #radius of the pellet
H = 1.0 #height of the pellet
dr = 1.1
n = int(R/dr)  # The pellet is divided in n annuli
print('n',n)
q = 450.0  # LHR of pin (kW/m)
Nc = q/2.0  # No of radial cracks
S = np.ones(n+1)  # initialization suRf_dotace area
r = np.ones(n+1)
r[0] = R
j = 2.0
V = float((3.14*(R**2)*H)/n)
print('V',V)

# Program

for i in range(n):
    r[i+1] = float(r[i]**2-((R**2)/n))**(0.5)
    print('i',i)
    S[0] = (2 * 3.14 * r[0] * H + 2 * 3.14 * (r[0]**2-r[1]**2)
            + 2 * H * (r[0] - r[1]) * Nc)
    if i<=j-1:
        S[i] = (2 * 3.14 * ((r[i])**2 - (r[i+1])**2)
            + 2 * H * (r[i]-r[i+1]) * Nc)
    elif i==j:
        S[i] = (2*3.14*r[i+1]*H + 2 * 3.14 * (r[i]**2 - r[i+1]**2)
            + 2 * H * (r[i]-r[i+1]) * Nc)
    elif i==j+1:
        S[i] = (2 * 3.14 * r[i+1] * H
                + 2 * 3.14 * (r[i+1]**2-r[i+2]**2))
    elif i<=n:
        S[i] = (2 * 3.14 * (r[i]**2 - r[i+1]**2))
print('r',r)
print('S',S)
print('S/V',S/V)
s_v_c = S/V
'''
a = np.array([1,2,3,5])
print a[-1]