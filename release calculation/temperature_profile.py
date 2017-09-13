'''
rf = pallet radious (meter)
dr = space between two grid point
alpha = 1/rho*Cp
dt = spacing between two time points (seconds)
sp = no of grid point
time _step = no of time iteration
kc = clad conductance (W/m-k)
G = gap between clad and fuel (meter)
Q  = heat generation term
tc = clad thickness (meter)
T = temperature in fuel pin (K)

'''

def temp_profile():
    import numpy as np
    from matplotlib import pyplot as plt
    import math
    from inputs_for_sv_cracks import input_data

    T = [397-100 for i in range(input_data.sp + 1)]

    # fuel surface temperature
    T[-1] = 330.0 + 273.0

    for j in range(0, input_data.time_step):
        for i in range(1, input_data.sp):
            T[1] = T[2]
            input_data.t = T[i] / 1000.0
            input_data.one = (1.0 / (input_data.A + input_data.B * input_data.t))
            input_data.two = (6400.0 / ((input_data.t)**(2.5)))
            input_data.threev = math.exp(-16.35 / (input_data.t))
            input_data.k[i] = 1.158*((1.0/float(input_data.A+input_data.B*input_data.t)) + (6400.0/(float(input_data.t)**(2.5))* math.exp(-16.35/float(input_data.t))))

            T[i] = (T[i] +
                    ((input_data.alphat * input_data.dt) / (i * (input_data.dr ** 2)))
                    * (((2 * i + 1) / 2) * input_data.k[i+1] * (T[i + 1] - T[i]) - ((2 * i - 1) / 2)
                       * input_data.k[i-1] * (T[i] - T[i - 1]))
                    + (input_data.Qe * input_data.alphat * input_data.dt))

    # print('T[-1]', T[-1])

    # Tc1 = T[-1] - (Qe * rf) / (hg * 2)
    # Tc2 = Tc1 - ((Qe * rf * tc) / (kc * 2))
    # pallet_temp = [T[-1, -1], Tc1, Tc2]
    # #print('pallet temp ', pallet_temp)

    ##print(k)

#    plt.plot(input_data.space, T[1:input_data.sp])
#    plt.show()
    return (T)
T = temp_profile()
# print(T)
#