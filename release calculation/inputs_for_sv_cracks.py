class input_data:
    ## Value Initilization
    import numpy as np
    from matplotlib import pyplot as plt
    import math



    # input data for temperature profile calculation

    rf = (5.52/2)*1E-3
    dr = 1E-4  # mesh spacing
    alphat = 1E-6  # alpha only contains 1/rho*Cp
    dt_max = (dr ** 2) / (2 * alphat)  # stability condition
    dt = dt_max / 10  # time step size
    sp = int(rf/dr) # no of spatial points
    time_step = int(1E4)  # time points
    kc = 29.0  # clad conductance
    kg = 58.22E-3  # conductance in gap between clad and gap
    G = 427E-6  # gap
    hg = kg / G  # gap conductance
    linear_heat_rate = 45E3
    Qe = linear_heat_rate/(3.14*rf**2)
    # Parameter for calculation of k
    k = [2 for i in range(sp + 1)]  # initial array of k
    x = 2.0 - 1.97  # 1.97 is O/M ratio
    A = 2.85 * x + 0.035
    B = -0.715 * x + 0.286
    space = [i for i in range(1, sp)]





    # Input data for crack surface area calculation

    a = 5E-6  # grain radius
#    R = rf-0.09E-3  # radius of the pellet
    R = (5.52/2)*1E-3
    H = 14.0E-3  # height of the pellet
    thickness_of_annuli = 0.5E-3
    no_of_annuli = int(R / thickness_of_annuli)  # The pellet is divided in n annuli
    q = 45 # E3 # LHR of pin (kW/m)
    Nc = q / 2.0  # No of radial cracks
    S = np.zeros(no_of_annuli + 1)  # initialization suRf_dotace area
    r = np.ones(no_of_annuli + 1)
    r[0] = R
    Vol = float((3.14 * (R ** 2) * H) / no_of_annuli)
    temp_T = np.ones(no_of_annuli)
    n = int(5E4)
    s_v_gb_cummulative = 0
    s_v_gb = np.zeros(n)
    fc = np.ones(n)
    Dfc = np.ones(n)
    Q = np.ones(n)
    DQ = np.ones(n)
    Rf_dot = np.ones(n)



    # time required for the establishment of grain bubble t_est

    # n = int(1E2)
    fc = np.ones(n)
    Dfc = np.ones(n)
    DQ = np.ones(n)
    Rf_dot = np.zeros(n)
    to = 1 * np.zeros(n)
    tc = 1E2*np.ones(n)
    rt = 1E-10*np.ones(n)
    # rt[1] = 1E-/7
    rt[1] = 1E-09
    rc = np.ones(n)
    rc[1] = 1E-7
    Re_dot = np.zeros(n)
    Re_dot[1] = 1E1
    Rc_dot = np.zeros(n)
    Rc_dot[1] = 1E1
    # rt[1] = 1E-7
    # rc[1] = 1E-7
    a1 = 0.1
    a2 = 2.2
    omega = 4.1E-29
    Del_S_by_S = np.zeros(n)
    N_dot = np.zeros(n)
    N = np.zeros(n)
    phi1 = np.zeros(n)
    u1 = np.zeros(n)
    Rc_dot[1] = 1
    tc_total = 1E1*np.ones(n)
    # N_max = np.ones(n)
    # N_max = np.zeros(n)


    # data for grain area/volume calculation

    alpha = 1E15   # fixed sink length (m-2)
    s = 3E-10  # atomic jump distance
    R_gas_constant = 1.987 # 1.9872036E-3 kcal/kelvin.mol  # Gas constant
    Z = 100  # n of sites around  a point from which recombination is inevitable
    rf1 = 0.5E-6
    K = 2E-4  # defect production rate per atom
    F = 1E15  # fission rate
    gamma = 0.5  # (check the value) #free suRf_dotace energy
    theta = 50
    fb = 0.25
    beta = 1.65E18  # (check the value)generation rate of stable fission gas atoms(value checked..!!)
    Pext = 0  # (check the value)  external pressure on the fuel pellet (value checked..!!)
    kb = 1.3807E-23  # J/kalvin#5.67E8 check value of boltzman constant
    bvd = 1E-13  # bvd = bv*delta


    # inputs for sv_cracks

#    lmbda = np.array ([2.05E-9, 6.7401E-7, 1.5424E-6, 3.6624E-6, 2.12E-5, 4.415E-5, 6.73E-5, 1.05E-4, 1.515E-4,
#                              7.35E-4, 7.55E-4, 8.17E-4, 3.01E-3, 3.75E-3])


#    lmbda = np.array([9.999999999999999547e-07,1.529999999999999963e-06,
#                      9.259999999999999433e-06,2.120000000000000031e-05,
#                      2.909999999999999921e-05,4.300000000000000165e-05,
#                      6.779999999999999460e-05,8.370000000000000161e-05,
#                      1.520000000000000084e-04,2.200000000000000078e-04,
#                      7.560000000000000530e-04,8.180000000000000422e-04,
#                      3.020000000000000115e-03])

    time = np.array([86400, 8280, 74880, 3150, 23652, 65188800, 948287520, 10668, 3392064, 3.22e7, 4365792, 907921440, 337435200, 16128, 4578, 10224, 452736, 32904, 583200, 7.6e11, 2.06e11, 451595520])

    lmbda = 0.693/time





    isotopes = ['I131', 'I132', 'I133', 'I134', 'I135', 'Cs134', 'Cs137', 'Rb88', 'Ru103', 'Ru-106', 'Sr-89', 'Sr90', 'Kr85', 'Kr85m', 'Kr87', 'Kr88', 'xe133', 'xe135', 'u237', 'Pu239', 'Pu240', 'Pu241']

#    isotopes = ['I-131', 'xe-133', 'I-133', 'xe-135', 'I-135', 'kr-85m',
#                'kr-88', 'I-132', 'kr-88', 'I-134',
#                'xe-135m', 'xe-138', 'xe-137']



    Ef = 3.2E-11
    L = 70E-3
    q = 40
    Y = np.array([0.0238,0.0670,0.0670,0.0654,0.0630,
     0.0130,0.0355,0.0430,0.0252,0.0773,0.0110,0.0642,0.0606])

    correction = np.array([1.000,1.205,1.100,3.353,1.000,1.019,1.002,
     11.70,1.016,2.100,6.373,1.003,1.072])

    production = (L*q*Y*correction)/Ef

    grainarea = np.ones(temp_T.size)
    Diff = np.ones(temp_T.size)
    rba = np.zeros(temp_T.size)
    s_v_t = np.zeros(temp_T.size)
    rbf = np.zeros(lmbda.size)



    # Information of isotopes
    no_of_iteration = int (1E6)
    dt_nu = 5E-1
    pin_failure_time = int (5E10)
    reactor_shut_down_time = 2E10




#    isotopes = ['kr-85', 'xe-131m', 'xe-133', 'xe-133m', 'xe-135', 'kr-85m', 'kr-88', 'kr-83m', 'kr-87', 'xe-134m',
#                'xe-135m', 'xe-138', 'xe-137', 'kr-89']

#    production =   1E12 * np.array([0.05, 0.09, 17.30, 0.57, 18.70, 1.40, 3.20, 0.74, 2.50, 0.47, 4.45, 12.80, 14.80, 3.7])

    # production_matrix = np.zeros([lmbda.size,lmbda.size])
    # lmbda_matrix = np.zeros([lmbda.size,lmbda.size])
    sweep_to_volume =   (60*3000*1E-6)/(4) # 2000 * 60 conversion of 2 litre/ min to cc/hout
    ##(200E-6 / (4))  # sweep rate = 0.5 m3/s, volume of the covergas=5 m3  (200E-6)




    Na = np.zeros([lmbda.size])

    Nb = np.zeros([lmbda.size])

    Na_inf = np.ones(lmbda.size)

    Nb_inf = np.ones(lmbda.size)

    # nu = 100





    # input date for nu calculations


    rupture_dia = 5E-3 # in m
    pin_dia = 5.1E-3

    area_of_rupture = 3.14 * (rupture_dia ** 2 / 4)  # Area of rupture
    height_plenum = 50E-3 # actual 320mm
    g = 9.8
    cde = 0.32
    gmma = 5.0 / 3
    c_gamma = (gmma * (2 / (gmma + 1)) ** ((gmma + 1) / (gmma - 1))) ** 0.5
    R_nu_gas = 08.314
    mu = 23.2
    plenum_volume = 3.14 * ((pin_dia+0.45E-3) ** 2 / 4) * height_plenum
    length_of_clad = 0.45E-3
    Tg = 400
    P1=0
    # P = 2E6*np.ones([no_of_iteration])
    # P1 = np.zeros([lmbda.size])
    # L = np.zeros([lmbda.size])
    # constant = np.zeros([lmbda.size])
    P_ref = 100
    L = 0
    isotope_index = 1
    '''
    ['kr-85'    0
    'xe-131m'   1
    'xe-133'    2
    'xe-133m'   3
    'xe-135'    4
    'kr-85m'    5
    'kr-88'     6
    'kr-83m'    7
    'kr-87'     8
    'xe-134m'   9
    'xe-135m'   10
    'xe-138'    11
    'xe-137'    12
    'kr-89'     13
    '''
