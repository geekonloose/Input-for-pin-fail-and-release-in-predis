def ugr_area(T):
    import numpy as np
    from matplotlib import pyplot as plt
    import math
    from inputs_for_sv_cracks import input_data

    # # print('T',T)

    # input_data.jv = 1E13 * math.exp(-5.42E4 / (input_data.R_gas_constant * T))  # vacancy jump rate
    input_data.jv = 7.6E-10 * math.exp(-7E4 / (1.986 * T))
    # # print('input_data.jv',input_data.jv)
    input_data.V1 = (input_data.alpha * (input_data.s ** 2) / (2 * input_data.Z)) * (
        ((1 + (4 * input_data.K * input_data.Z) / (input_data.jv * input_data.alpha * (input_data.s ** 2))) ** 0.5) - 1)
    # # print('input_data.V1',input_data.V1)
    D = 7.6E-10 * math.exp(-7E4 / (input_data.R_gas_constant * T)) + (
                                                                     input_data.s ** 2) * input_data.jv * input_data.V1 + 2E-40 * input_data.F
    # D = D
    # # print("D", D)

    input_data.Ds = float(1.12E3 * math.exp(-119000 / (input_data.R_gas_constant * T)))
    input_data.f_theta = 1 - (1.5) * math.cos(input_data.theta) + 0.5 * (math.cos(input_data.theta)) ** 3
    input_data.N_max = ((4 * input_data.rf1 * input_data.f_theta * input_data.fb) / (
        3 * input_data.kb * T * (math.sin(input_data.theta)) ** 2)) * (
                           ((2 * input_data.gamma) / (
                               input_data.rf1)) + input_data.Pext)  # area density of fission gas atoms on grain faces at saturation

    # calculation of delay at grain face bubble




    i = 1
    while i < input_data.n:

        # # print('i', i)

        input_data.rt[i] = input_data.rt[i-1]*1.12*input_data.a2
        # print 'inputy_data.rt',input_data.rt[i]
        # input_data.to[i] = input_data.a1*(input_data.kb*T/
        #                                   (input_data.Ds*input_data.gamma*(input_data.omega)**(4.0/3)))*input_data.rt[i]**4
        input_data.to[i] = (input_data.a1 * (
            (input_data.kb * T) / (input_data.Ds * input_data.gamma * (input_data.omega ** (float(4 / 3))))) * (
                                           input_data.rt[i] ** 4))



        # for determining the rf_dot
        input_data.w = (D * i) / (input_data.a ** 2)

        if input_data.w <= (3.14 ** (-2)) and i > 0:

            input_data.fc[i] = 4 * (input_data.w / 3.14) ** 0.5 - 1.5 * input_data.w
            input_data.Dfc[i] = 2 * (D / ((input_data.a ** 2) * i * 3.14)) ** 0.5 - 1.5 * (
                D / (input_data.a ** 2))

        elif input_data.w >= (3.14 ** (-2)) and i > 0:

            input_data.fc[i] = 1 - (0.0662 / input_data.w) * (
                1 - 0.93 * math.exp(
                    (-3.14 ** (-2)) * input_data.w))  # for w<=3.14**-2#fractional release of stable fission gas atoms
            input_data.Dfc[i] = ((0.662 * (input_data.a ** 2) / (D * i ** 2))
                                            - (0.93 * 0.0662 * (input_data.a ** 2) / (
                D * (i ** 2)) + 0.93 * 0.0662 / ((3.14 ** 2) * i))
                                            * math.exp(-D * i / ((3.14 ** 2) * input_data.a ** 2)))

        input_data.DQ[i] = input_data.fc[i] + i * input_data.Dfc[i]


        input_data.Rf_dot[i] = (2.0 * input_data.a / 3.0) * input_data.DQ[i]

        # calculation of the N[i]


        input_data.u1[i] = ((4 * input_data.bvd * (i) ** (0.5) / (D * 3.14) ** 0.5))

        input_data.phi1[i] = input_data.A * (
            input_data.u1[i] ** 2 - 2 * input_data.u1[i] + 2 - 2 * math.exp(
                -input_data.u1[i]))

        input_data.N[i] = input_data.phi1[i] * input_data.N_max


        input_data.N_dot[i] = (4 * input_data.beta * (D * i / 3.14) ** (0.5) * (
            1 - ((input_data.bvd) * input_data.N[i]) / (2 * D * input_data.beta * i)))



        # # print(input_data.rc[i])
        # loop
        i3 = 1
        while 2==1+1:

            input_data.Del_S_by_S[i3] = input_data.Del_S_by_S[i3-1]+0.541*(input_data.rc[i3]/input_data.a)
            # # print(input_data.Del_S_by_S[i3])
            input_data.Rc_dot[i3] = input_data.Rc_dot[i3-1]+1.33*input_data.rc[i3]*(input_data.Del_S_by_S[i3]*input_data.N_dot[i3]
                                                            + input_data.Rf_dot[i3])
            # # print('input_data.Rc_dot[i3]',input_data.Rc_dot[i3])
            input_data.Re_dot[i3] = input_data.Re_dot[i3-1]+0.834*input_data.a*(input_data.Del_S_by_S[i3]*input_data.N_dot[i3])


            # # print('input_data.Re_dot[i3]',input_data.Re_dot[i3])
            # input_data.tc[i3] = input_data.tc[i3-1]+(0.115*3.14)/(input_data.kb*T)\
            #                    *((input_data.a2*input_data.rt[i3])**2/input_data.Re_dot[i3])\
            #                    *(1.34*input_data.gamma/(input_data.a2*input_data.rt[i3]))
            # # print('input_data.tc[i3]',input_data.tc[i3])
            # # print(i)
            input_data.N_max = (((0.0194 * 3.14 * (input_data.a) ** 3) / (input_data.kb * T)) *
                                   (((2.27 * input_data.gamma )/ input_data.a) ))
            # # print(input_data.Rc_dot[i])
            # # print(input_data.N_max)
            input_data.tc_total[i3] = input_data.tc_total[i3-1]+(input_data.N_max)/(input_data.Rc_dot[i3])

            # # print 'input_data.tc_total',input_data.tc_total[i3]

            # print i3
            input_data.rc[i3] = 0.4965 * input_data.a + (1.54 * input_data.a * i3) / (
            3 * input_data.to[i3] + 4 * input_data.tc[i3])
            # # print('input_data.to[i3]',input_data.tc[i3])
            # # print 'input_data.rc[i3]', input_data.rc[i3]
            if (input_data.rc[i3] - input_data.rc[i3 - 1] <= 1E-10):
                break

            i3=i3+1
        # input_data.rc1[i] = input_data.rc[i3]
        # print 'rt',input_data.rt[i]
        if input_data.rt[i]>(0.72*input_data.a)/(2*input_data.a2):
            break
        i = i+1

    rc1 = input_data.rc[i3]
    input_data.s_v_gb = 1.623*(rc1**2)/(input_data.a**3)
    return (input_data.s_v_gb,D)
    print('gb',s_v_gb)
    # # print('rc',input_data.rc[i3])