
def nu_modelling(j1):
    import numpy as np
    from matplotlib import pyplot as plt
    import math
    from inputs_for_sv_cracks import input_data
    # j1=100

    # isotopes = ['kr-85', 'xe-131m', 'xe-133', 'xe-133m', 'xe-135', 'kr-85m', 'kr-88', 'kr-83m', 'kr-87', 'xe-134m',
    #             'xe-135m', 'xe-138', 'xe-137', 'kr-89']
    # lmbda = np.array([2.05E-9, 6.7401E-7, 1.5424E-6, 3.6624E-6, 2.12E-5, 4.415E-5, 6.73E-5, 1.05E-4, 1.515E-4, 7.35E-4, 7.55E-4, 8.17E-4,
    #      3.01E-3, 3.75E-3])

    # for i1 in range(lmbda.size-1):
    # for j1 in range(no_of_iteration-1):
    input_data.t_nu[j1] = j1*input_data.dt_nu
    input_data.constant[j1] = (input_data.rupture_dia**2*input_data.g*input_data.area_of_rupture*input_data.rho[j1])/(64*input_data.length_of_clad*input_data.mu)
    input_data.P[j1+1] = input_data.P[j1] - input_data.constant[j1]*(input_data.P[j1]-input_data.P_ref)*input_data.P[j1]*input_data.dt_nu
    input_data.rho[j1+1] = input_data.P[j1+1]/(input_data.R_nu_gas*input_data.Tg)
    # input_data.pc[j1+1]input_data.gmma




    input_data.L[j1] = (input_data.constant[j1]*(input_data.P[j1]-input_data.P_ref))

    input_data.L[j1] = (1/input_data.plenum_volume)*input_data.L[j1]
    input_data.L[j1] = (1/input_data.rho[j1])*input_data.L[j1]
    return input_data.L[j1]


