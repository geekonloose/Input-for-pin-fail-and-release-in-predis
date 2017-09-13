from inputs_for_sv_cracks import input_data


def gr_area(T):
    import numpy as np
    from matplotlib import pyplot as plt
    import math
    from inputs_for_sv_cracks import input_data
    # print('T',T)
    
    # input_data.jv = 1E13 * math.exp(-5.42E4 / (input_data.R_gas_constant * T))  # vacancy jump rate
    input_data.jv = 7.6E-10* math.exp(-7E4/(1.986*T))
    print('input_data.jv',input_data.jv)
    input_data.V1 = (input_data.alpha * (input_data.s ** 2) / (2 * input_data.Z)) * (
    ((1 + (4 * input_data.K * input_data.Z) / (input_data.jv * input_data.alpha * (input_data.s ** 2))) ** 0.5) - 1)
    # print('input_data.V1',input_data.V1)
    D = 7.6E-10 * math.exp(-7E4 /(input_data.R_gas_constant* T)) + (input_data.s ** 2) * input_data.jv * input_data.V1 + 2E-40 * input_data.F
    # D = D
    # print("D", D)

    input_data.Ds = float(1.12E3 * math.exp(-119000 / (input_data.R_gas_constant * T)))
    input_data.f_theta = 1 - (1.5) * math.cos(input_data.theta) + 0.5 * (math.cos(input_data.theta)) ** 3
    input_data.N_max = ((4 * input_data.rf1 * input_data.f_theta * input_data.fb) / (
    3 * input_data.kb * T * (math.sin(input_data.theta)) ** 2)) * (
                           ((2 * input_data.gamma) / (
                           input_data.rf1)) + input_data.Pext)  # area density of fission gas atoms on grain faces at saturation

    # calculation of delay at grain face bubble




    input_data.i = 1
    while input_data.i < input_data.n:

        if input_data.rt[input_data.i] < (0.75 * input_data.a) / 2 * input_data.a2:
            input_data.rt[input_data.i] = 1.12 * input_data.a2 * input_data.rt[input_data.i - 1]
            input_data.to[input_data.i] = (input_data.a1 * (
            (input_data.kb * T) / (input_data.Ds * input_data.gamma * (input_data.omega ** (float(4 / 3))))) * (
                                           input_data.rt[input_data.i] ** 4))
        elif input_data.rt[input_data.i] > (0.75 * input_data.a) / 2 * input_data.a2:
            input_data.rt[input_data.i] = 0.466 * ((input_data.rc[input_data.i] ** 3) / input_data.a) ** 0.5
            input_data.to[input_data.i] = (input_data.a1 * input_data.kb * T * (0.75 * input_data.a) ** 2) / (
            input_data.Ds * input_data.gamma * (input_data.omega) ** (4.0 / 3) * (2 * input_data.a2) ** 2) * \
                                          input_data.rt[input_data.i] ** 2

        if input_data.i < 98:
            input_data.rc[input_data.i] = 0.4965 * input_data.a + (1.54 * input_data.i) / (
            3 * input_data.to[input_data.i] + 7 * input_data.tc_total)


        # calculate del_s_by_s
        input_data.Del_S_by_S[input_data.i] = 0.541 * (input_data.rc[input_data.i] / input_data.a)

        # for determining the rf_dot
        input_data.w = (D * input_data.i) / (input_data.a ** 2)

        if input_data.w <= (3.14 ** (-2)) and input_data.i > 0:

            input_data.fc[input_data.i] = 4 * (input_data.w / 3.14) ** 0.5 - 1.5 * input_data.w
            input_data.Dfc[input_data.i] = 2 * (D / ((input_data.a ** 2) * input_data.i * 3.14)) ** 0.5 - 1.5 * (
            D / (input_data.a ** 2))

        elif input_data.w >= (3.14 ** (-2)) and input_data.i > 0:

            input_data.fc[input_data.i] = 1 - (0.0662 / input_data.w) * (
                1 - 0.93 * math.exp(
                    (-3.14 ** (-2)) * input_data.w))  # for w<=3.14**-2#fractional release of stable fission gas atoms
            input_data.Dfc[input_data.i] = ((0.662 * (input_data.a ** 2) / (D * input_data.i ** 2))
                                            - (0.93 * 0.0662 * (input_data.a ** 2) / (
            D * (input_data.i ** 2)) + 0.93 * 0.0662 / ((3.14 ** 2) * input_data.i))
                                            * math.exp(-D * input_data.i / ((3.14 ** 2) * input_data.a ** 2)))

        input_data.DQ[input_data.i] = input_data.fc[input_data.i] + input_data.i * input_data.Dfc[input_data.i]

        # calculation of Rf_dot
        input_data.Rf_dot[input_data.i] = (2.0 * input_data.a / 3.0) * input_data.DQ[input_data.i]



        # calculation of the N[input_data.i]


        input_data.u1[input_data.i] = ((4 * input_data.bvd * (input_data.i) ** (0.5) / (D * 3.14) ** 0.5))

        input_data.phi1[input_data.i] = input_data.A * (
        input_data.u1[input_data.i] ** 2 - 2 * input_data.u1[input_data.i] + 2 - 2 * math.exp(
            -input_data.u1[input_data.i]))

        input_data.N[input_data.i] = input_data.phi1[input_data.i] * input_data.N_max

        # calculation of N_dot

        input_data.N_dot[input_data.i] = (4 * input_data.beta * (D * input_data.i / 3.14) ** (0.5) * (
        1 - ((input_data.bvd) * input_data.N[input_data.i]) / (2 * D * input_data.beta * input_data.i)))


        # find Rc_dot
        input_data.Rc_dot[input_data.i] = 1.133 * (input_data.rc[input_data.i] ** 2) * (
        input_data.Del_S_by_S[input_data.i] * input_data.N_dot[input_data.i] + input_data.Rf_dot[input_data.i])


        # find Nc_max
        input_data.Nc_max = ((0.0194 * 3.14 * (input_data.a) ** 3) / (input_data.kb * T) * (
        (2.27 * input_data.gamma / input_data.a) + input_data.Pext))



        # find tc_total
        input_data.tc_total = input_data.Nc_max / input_data.Rc_dot[input_data.i]

        if input_data.to[input_data.i] > 5E8:
            break


        input_data.i = input_data.i + 1

    for num in range(input_data.n):
        if input_data.rt[num] == 1:
            input_data.rt[num] = 0
            input_data.rc[num] = 0

    input_data.test = sum(input_data.to[1:input_data.i]) + sum(input_data.tc[1:input_data.i])

    input_data.p = np.array([num for num in range(input_data.rc.size)])
    input_data.rc = input_data.rc / input_data.a
    print(input_data.rc)
    input_data.s_v_gb = (input_data.Del_S_by_S[input_data.i] * 4 * 3.14 * input_data.a ** 2) / (
    (4.0 / 3) * 3.14 * input_data.a ** 3)


    return (input_data.s_v_gb, D)
gr_area(2000)
