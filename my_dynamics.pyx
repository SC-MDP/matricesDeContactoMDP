import random

def evol_sus(double p, int stop):
    cdef int sus = 0
    cdef int exp = 0
    cdef double ran = 0.0
    cdef int i = 0

    if p < 0.000000000001:
        i = stop
    
    while i < stop:
        i += 1
        ran = random.random()
        if ran < p:
            exp += 1
    sus -= exp

    return exp, sus

def evol_exp(double p_a, double p_e, int stop):
    cdef int exp = 0
    cdef int inf = 0
    cdef int asym = 0
    cdef double ran_a = 0.0
    cdef double ran_e = 0.0
    cdef int i = 0

    while i < stop:
        i += 1
        ran_a = random.random()
        ran_e = random.random()
        if ran_a < p_a:
            if ran_e < p_e:
                inf += 1
            else:
                asym += 1
    exp -= inf + asym

    return inf, asym, exp

def evol_inf(double p_o, double p_x, double t, double t_s, double t_c,\
             double t_cs, int stop, int flag_beds, int beds):
    cdef int inf = 0
    cdef int inf_acc = 0
    cdef int hosp = 0
    cdef int hosp_s = 0
    cdef int hosp_c = 0
    cdef int hosp_cs = 0
    cdef int dead = 0
    cdef int mild = 0
    cdef double ran_o = 0.0
    cdef double ran_x = 0.0
    cdef int i = 0

    while i < stop:
        i += 1
        ran_o = random.random()
        ran_x = random.random()
        if ran_o < p_o:
            inf_acc += 1
            if ran_x < p_x:             
                h_index = random.choices(['H','Hs','Hc','Hcs'],
                                        [1.0 - (t_s + t_c + t_cs), t_s, t_c, t_cs])
                if h_index[0]=='Hs' or h_index[0]=='Hcs':
                    if flag_beds >= beds:
                        dead += 1
                    else:
                        flag_beds += 1
                        if h_index[0] == 'Hs':
                            hosp_s += 1
                        else:
                            hosp_cs += 1
                else:
                    if h_index[0] == 'H':
                        hosp += 1
                    elif h_index[0] == 'Hc':
                        hosp_c += 1
            else:
                mild += 1
    inf -= inf_acc

    return inf, inf_acc, hosp, hosp_s, hosp_c,\
           hosp_cs, mild, dead, flag_beds

def evol_asym(double p_t, double p_g, int stop):
    cdef int asym = 0
    cdef int mild = 0
    cdef inf_acc = 0
    cdef int a_index = 0
    cdef int i = 0

    while i < stop:
        i += 1
        as_index = random.choices(['M','R','A'], [p_t, p_g, 1.0 - (p_t + p_g)])
        if as_index[0] == 'M':
            asym -= 1
            mild += 1
        elif as_index[0] == 'R':
            asym -= 1
        else:
            pass
    inf_acc = mild

    return mild, asym, inf_acc

def evol_hosp(double p, int stop):
    cdef int v_in = 0
    cdef int v_out = 0
    cdef double ran = 0
    cdef int i = 0

    while i < stop:
        i += 1
        ran = random.random()
        if ran < p:
            v_in += 1
    v_out -= v_in

    return v_in, v_out

def evol_mild(double p, int stop):
    cdef int v_in = 0
    cdef int v_out = 0
    cdef double ran = 0
    cdef int i = 0

    while i < stop:
        i += 1
        ran = random.random()
        if ran < p:
            v_in += 1
    v_out -= v_in

    return v_in, v_out
    

