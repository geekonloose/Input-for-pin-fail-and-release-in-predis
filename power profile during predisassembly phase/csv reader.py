# -*- coding: utf-8 -*-
"""
Created on Sat Jul  1 19:18:34 2017

@author: parth

Test file for reading csv files
"""

import csv
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import math
x = []
y = []
with open('power_profile.csv') as csvfile:
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
print(x1)
print(y)
def func(x1,a,b,c,d,e,f):
    return a* x1**5 + b*x1**4 + c*x1**3 + d*x1**2+e*x1+f


params = curve_fit(func, x1, y)
[a, b,c,d,e,f] = params[0]
plt.plot(x1,y,color='red')

y = a* x1**5 + b*x1**4 + c*x1**3 + d*x1**2+e*x1+f
plt.plot(x1,y,color='blue')

plt.show()
