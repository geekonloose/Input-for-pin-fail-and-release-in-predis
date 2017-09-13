# Pre-prossessing of temperature profile

from inputs_for_sv_cracks import input_data
def crack():
    from inputs_for_sv_cracks import input_data

    input_data.j = input_data.no_of_annuli-1

    for i in range(input_data.no_of_annuli):
        input_data.r[i + 1] = float(input_data.r[i] ** 2 - ((input_data.R ** 2) / input_data.no_of_annuli)) ** 0.5
        # #print('i', i)
        input_data.S[0] = (
        2 * 3.14 * input_data.r[0] * input_data.H + 2 * 3.14 * (input_data.r[0] ** 2 - input_data.r[1] ** 2)
        + 2 * input_data.H * (input_data.r[0] - input_data.r[1]) * input_data.Nc)
        if i <= input_data.j - 1:
            input_data.S[i] = (2 * 3.14 * ((input_data.r[i]) ** 2 - (input_data.r[i + 1]) ** 2)
                               + 2 * input_data.H * (input_data.r[i] - input_data.r[i + 1]) * input_data.Nc)
        elif i == input_data.j:
            input_data.S[i] = (
            2 * 3.14 * input_data.r[i + 1] * input_data.H + 2 * 3.14 * (input_data.r[i] ** 2 - input_data.r[i + 1] ** 2)
            + 2 * input_data.H * (input_data.r[i] - input_data.r[i + 1]) * input_data.Nc)
        elif i == input_data.j + 1:
            input_data.S[i] = (2 * 3.14 * input_data.r[i + 1] * input_data.H
                               + 2 * 3.14 * (input_data.r[i + 1] ** 2 - input_data.r[i + 2] ** 2))
        elif i <= input_data.no_of_annuli:
            print(i)
            print('ri',input_data.r[i])
            print('ri+1',input_data.r[i+1])
            input_data.S[i] = (2 * 3.14 * (input_data.r[i] ** 2 - input_data.r[i + 1] ** 2))

    input_data.s_v_c = input_data.S / input_data.Vol

    return ()


crack()
print(input_data.s_v_c)
