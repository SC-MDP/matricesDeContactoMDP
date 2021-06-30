import my_dynamics as md
import math
import random
import numpy as np
import sys
import os
##sys.exit("STOPPED BY CODE")

#Clase grupo etario.
class Grupo:
    pass

#Parametros del programa y variables auxiliares
nrea = 3
batch = 10
n_days = 31
dt = 1/batch
time = 0
q0 = 0.5
tass = 7
beds  = 60
r0 = 0.02
aux_r0 = r0
eps_a = [0.333 if i < 4 else 0.333 if 4 <= i < 18 else 1 for i in range(21)]
#Lectura parametros epidemicos por rango etario
zeta = np.loadtxt('ParamEtar/p_zeta.txt')
thet = np.loadtxt('ParamEtar/p_thet.txt')
thet_st = np.loadtxt('ParamEtar/p_thet_st.txt')
thet_cr = np.loadtxt('ParamEtar/p_thet_cr.txt')
thet_cr_st = np.loadtxt('ParamEtar/p_thet_cr_st.txt')

#Lectura semillas iniciales por compartimento
input_ajuste = np.loadtxt('ParamEtar/input.txt',usecols=range(1,22),max_rows=11)

S_input = input_ajuste[0]
E_input = input_ajuste[1]
A_input = input_ajuste[2]
I_input = input_ajuste[3]
M_input = input_ajuste[4]
HnL_input = input_ajuste[5]
HsL_input = input_ajuste[6]
HnD_input = input_ajuste[7]
HsD_input = input_ajuste[8]
R_input = input_ajuste[9]
D_input = input_ajuste[10]
Iacc_input = np.loadtxt('ParamEtar/input.txt',usecols=range(1,2),skiprows=11)

#Lectura de las matrices de contacto
c_work = np.loadtxt('Matrices/c_trabajo.dat')
c_educ = np.loadtxt('Matrices/c_educacion.dat')
c_hous = np.loadtxt('Matrices/c_hogares.dat')
c_rand = np.loadtxt('Matrices/c_aleatoria.dat')
#Lectura de datos poblacionales
pob = np.loadtxt('Poblacion/Pob5_25.txt', dtype=np.int32)

grupos_etar = c_work.shape[0] #Numero de grupos etarios.

c_matrixOut = np.zeros((grupos_etar, grupos_etar)) #matriz de los esenciales
c_matrixIn = np.zeros((grupos_etar, grupos_etar)) #matriz del grupo de no-esenciales aislados 

c_list = [c_work, c_educ, c_hous, c_rand] #Lista de matrices
c_weights = [0.19, 0.18, 0.30, 0.33] #Normalizado. Peso de los esenciales y los que no se aislan temporalmente

for i in range(len(c_list)): #Suma pesada de las matrices
    c_matrixOut = np.add(c_matrixOut, c_weights[i] * c_list[i])

#Loop de realizaciones
for rea in range(nrea):
    print('rea', rea)
    
    #Lista de redes. Cada red es un grupo etario
    G = []
    Ntotal = 0
    flag_stop = 0
    flag_beds = 0
    Inf_init = 0.0

    for i in range(grupos_etar):
        G.append(Grupo())
        G[i].N = int(pob[i] * 1.060529347048)
        
        #Subdivisión del grupo. "0" son los esenciales, A y B son los grupos que se alternan
        G[i].NumA = int(round(G[i].N * (1 - q0) / 2, 0))
        G[i].NumB = int(round(G[i].N * (1 - q0) / 2, 0))
        G[i].Num0 = G[i].N - G[i].NumA - G[i].NumB
     
        #Parametros epidémicos del grupo
##        G[i].beta = bet_up[i] 
##        G[i].beta_a = bet_up[i]
        G[i].beta = 5.6
        G[i].beta_a = 5.6
        G[i].alpha = 1.0/(5.1-2.0)
        G[i].epsilon = eps_a[i]
        G[i].omega = 1.0/(4.0+2.0)
        G[i].xi = zeta[i]
        G[i].tita = thet[i]
        G[i].tita_star = thet_st[i]
        G[i].tita_cross = thet_cr[i]
        G[i].tita_cr_st = thet_cr_st[i]
        G[i].gamma = 1.0/8.0
        G[i].gamma_mild = 1.0/14.0
        G[i].gamma_asym = 1.0/9.5
        G[i].gamma_star = 1.0/10.0
        G[i].delta = 1.0/8.0
        G[i].delta_star = 1.0/11.0

        #Nro de casos segun estado de salud dentro del grupo        
        G[i].n_sus = int(round(S_input[i],0))
        G[i].n_exp = int(round(E_input[i],0))
        G[i].n_inf = int(round(I_input[i],0))
        G[i].n_asym = int(round(A_input[i],0))
        G[i].n_hosp =  int(round(HnL_input[i],0))
        G[i].n_hosp_star = int(round(HsL_input[i],0))
        G[i].n_hosp_cross = int(round(HnD_input[i],0))
        G[i].n_hosp_cr_st = int(round(HsD_input[i],0))
        G[i].n_mild = int(round(M_input[i],0))
        G[i].n_rec = int(round(R_input[i],0))
        G[i].n_dead = int(round(D_input[i],0))
        
        G[i].n_inf_acc = int(round(Iacc_input/grupos_etar,0)) #No importa distribucion porque sacamos la suma

        #Variable auxiliar
        Ntotal += G[i].N
        Inf_init += G[i].n_inf
        flag_beds += G[i].n_hosp_star + G[i].n_hosp_cr_st

    #Estados extra por estrategia
    for g in range(grupos_etar):  
        G[g].n_inf0 = G[g].n_infA = G[g].n_infB = 0        
        G[g].n_susA = int(round(G[g].n_sus * (1 - q0) / 2, 0))
        G[g].n_susB = int(round(G[g].n_sus * (1 - q0) / 2, 0))
        G[g].n_sus0 = G[g].n_sus - G[g].n_susA - G[g].n_susB
        G[g].n_exp0 = G[g].n_expA = G[g].n_expB = 0	
        G[g].n_asym0 = G[g].n_asymA = G[g].n_asymB = 0
        G[g].b_eff_tot0 = G[g].b_eff_totA = G[g].b_eff_totB = 0

        #Registro las semillas en los subgrupos correspondientes
        for i in range(G[g].n_inf):
            ran_sub = np.random.choice(['0', 'A', 'B'], p=[G[g].Num0/G[g].N, G[g].NumA/G[g].N, G[g].NumB/G[g].N])
            for i in ran_sub:
                if i == '0' and G[g].n_sus0 > 0:
                    G[g].n_inf0 +=1
##                    G[g].n_sus0 -=1
                elif i == 'A' and G[g].n_susA > 0:
                    G[g].n_infA +=1
##                    G[g].n_susA -=1
                elif i == 'B' and G[g].n_susB > 0:
                    G[g].n_infB +=1
##                    G[g].n_susB -=1
                else:
                    print('Warning: algo no esta bien con las semillas I')

        for i in range(G[g].n_asym):
            ran_sub = np.random.choice(['0', 'A', 'B'], p=[G[g].Num0/G[g].N, G[g].NumA/G[g].N, G[g].NumB/G[g].N])
            for i in ran_sub:
                if i == '0' and G[g].n_sus0 > 0:
                    G[g].n_asym0 +=1
##                    G[g].n_sus0 -=1
                elif i == 'A' and G[g].n_susA > 0:
                    G[g].n_asymA +=1
##                    G[g].n_susA -=1
                elif i == 'B' and G[g].n_susB > 0:
                    G[g].n_asymB +=1
##                    G[g].n_susB -=1
                else:
                    print('Warning: algo no esta bien con las semillas A')

        for i in range(G[g].n_exp):
            ran_sub = np.random.choice(['0', 'A', 'B'], p=[G[g].Num0/G[g].N, G[g].NumA/G[g].N, G[g].NumB/G[g].N])
            for i in ran_sub:
                if i == '0' and G[g].n_sus0 > 0:
                    G[g].n_exp0 +=1
##                    G[g].n_sus0 -=1
                elif i == 'A' and G[g].n_susA > 0:
                    G[g].n_expA +=1
##                    G[g].n_susA -=1
                elif i == 'B' and G[g].n_susB > 0:
                    G[g].n_expB +=1
##                    G[g].n_susB -=1
                else:
                    print('Warning: algo no esta bien con las semillas E')

    #Listas output
    ICUtot = 0
    Itot = 0
    Dtot = 0
    Sustot = 0
    Asitot = 0
    Acutot = 0

    #Creo archivo. Se fija los nombres ya existentes.
    #Creo archivo.
    ffnn2 = '_q0'+str(format(q0, '03g'))+'_tass'+str(tass)+'_r0'+str(format(aux_r0, '03g'))
    if not os.path.exists('results_est%s' %ffnn2):
        os.makedirs('results_est%s' %ffnn2)

    i = 0
    ffnn = str(format(i, '04d'))
    while os.path.exists('results_est%s/EvolInfectados_%s.dat' %(ffnn2, ffnn)):
        i += 1
        ffnn = str(format(i, '04d'))
    
    file=open('results_est%s/EvolInfectados_%s.dat' %(ffnn2, ffnn),'w')
    fileED=open('results_est%s/EvolDead_%s.dat' %(ffnn2, ffnn),'w')
    fileICU=open('results_est%s/EvolHICU_%s.dat' %(ffnn2, ffnn),'w')
    fileSus=open('results_est%s/EvolSus_%s.dat' %(ffnn2, ffnn),'w')
    fileAsi=open('results_est%s/EvolAsimp_%s.dat' %(ffnn2, ffnn),'w')
    fileAcu=open('results_est%s/EvolAcum_%s.dat' %(ffnn2, ffnn),'w')

    #DINAMICA
    time = 0
    while flag_stop < Ntotal and time < n_days*batch:
        flag_stop = 0

        #Cantidades acumuladas
        Aux=0
        for g in range(grupos_etar):
            Aux = Aux + G[g].n_inf
        Itot = Aux   
            
        Aux=0
        for g in range(grupos_etar):
            Aux = Aux + G[g].n_hosp_star + G[g].n_hosp_cr_st
        ICUtot = Aux

        Aux=0
        for g in range(grupos_etar):
            Aux=Aux + G[g].n_dead
        Dtot = Aux

        Aux=0
        for g in range(grupos_etar):
            Aux = Aux + G[g].n_sus
        Sustot = Aux
        
        Aux=0
        for g in range(grupos_etar):
            Aux = Aux + G[g].n_asym
        Asitot = Aux
        
        Aux=0
        for g in range(grupos_etar):
            Aux = Aux + G[g].n_inf_acc
        Acutot = Aux
        
        #Salvo en los archivos
        if time%batch==0:
            file.write(str(time*dt) + "\t" + str(Itot) + "\n") 
            fileED.write(str(time*dt) + "\t" + str(Dtot) + "\n") 
            fileICU.write(str(time*dt) + "\t" + str(ICUtot) + "\n") 
            fileSus.write(str(time*dt) + "\t" + str(Sustot) + "\n") 
            fileAsi.write(str(time*dt) + "\t" + str(Asitot) + "\n") 
            fileAcu.write(str(time*dt) + "\t" + str(Acutot) + "\n")

        #Estrategia de tests
        if Itot <= Inf_init: #corto estrategia de tests
            r0 = 0.0
        else:
            r0 = aux_r0
        
        #Calculo de las fuerzas de infección
        if math.sin(math.pi*(time*dt/tass))>=0:
            for g in range(grupos_etar):
                G[g].fout = G[g].fin = 0
                b_eff_s_in = b_eff_a_in = b_eff_s_out = b_eff_a_out = 1 #auxiliares
                    
                for i in range(grupos_etar):
                    b_eff_s_out *= (1-G[g].beta * c_matrixOut[g][i] * dt / G[i].N) ** (G[i].n_inf0+G[i].n_infA) #estocasticas
                    b_eff_a_out *= (1-G[g].beta_a * c_matrixOut[g][i] * dt / G[i].N) ** (G[i].n_asym0+G[i].n_asymA) #estocasticas
                    
                G[g].b_eff_tot0 = G[g].b_eff_totA = 1 - b_eff_s_out * b_eff_a_out
                G[g].b_eff_totB = 0
                
        else:
            for g in range(grupos_etar):
                G[g].fout = G[g].fin = 0
                b_eff_s_in = b_eff_a_in = b_eff_s_out = b_eff_a_out = 1 #auxiliares
                    
                for i in range(grupos_etar):
                    b_eff_s_out *= (1-G[g].beta * c_matrixOut[g][i] * dt / G[i].N) ** (G[i].n_inf0+G[i].n_infB) #estocasticas
                    b_eff_a_out *= (1-G[g].beta_a * c_matrixOut[g][i] * dt / G[i].N) ** (G[i].n_asym0+G[i].n_asymB) #estocasticas
                    
                G[g].b_eff_tot0 = G[g].b_eff_totB = 1 - b_eff_s_out * b_eff_a_out
                G[g].b_eff_totA = 0
          
        #Evolución de la enfermedad    
        for g in range(grupos_etar):
            
            #Evolución de susceptibles
            aux_exp0, aux_expA, aux_expB = G[g].n_exp0, G[g].n_expA, G[g].n_expB #Guardo nros viejos para la actualizacion de este compartimento

            ev_exp, ev_sus = md.evol_sus(G[g].b_eff_tot0, G[g].n_sus0)
            G[g].n_exp += ev_exp
            G[g].n_exp0 += ev_exp
            G[g].n_sus += ev_sus
            G[g].n_sus0 += ev_sus

            ev_exp, ev_sus = md.evol_sus(G[g].b_eff_totA, G[g].n_susA)
            G[g].n_exp += ev_exp
            G[g].n_expA += ev_exp
            G[g].n_sus += ev_sus
            G[g].n_susA += ev_sus

            ev_exp, ev_sus = md.evol_sus(G[g].b_eff_totB, G[g].n_susB)
            G[g].n_exp += ev_exp
            G[g].n_expB += ev_exp
            G[g].n_sus += ev_sus
            G[g].n_susB += ev_sus            
                    
            #Evolución de expuestos
            aux_inf0, aux_infA, aux_infB = G[g].n_inf0, G[g].n_infA, G[g].n_infB
            aux_asym0, aux_asymA, aux_asymB = G[g].n_asym0, G[g].n_asymA, G[g].n_asymB

            ev_inf, ev_asym, ev_exp = md.evol_exp(G[g].alpha * dt, G[g].epsilon, aux_exp0)
            G[g].n_inf += ev_inf
            G[g].n_inf0 += ev_inf
            G[g].n_asym += ev_asym
            G[g].n_asym0 += ev_asym
            G[g].n_exp += ev_exp
            G[g].n_exp0 += ev_exp

            ev_inf, ev_asym, ev_exp = md.evol_exp(G[g].alpha * dt, G[g].epsilon, aux_expA)
            G[g].n_inf += ev_inf
            G[g].n_infA += ev_inf
            G[g].n_asym += ev_asym
            G[g].n_asymA += ev_asym
            G[g].n_exp += ev_exp
            G[g].n_expA += ev_exp

            ev_inf, ev_asym, ev_exp = md.evol_exp(G[g].alpha * dt, G[g].epsilon, aux_expB)
            G[g].n_inf += ev_inf
            G[g].n_infB += ev_inf
            G[g].n_asym += ev_asym
            G[g].n_asymB += ev_asym
            G[g].n_exp += ev_exp
            G[g].n_expB += ev_exp        

            #Evolución de infectados
            aux_hosp, aux_hosp_st, aux_hosp_cr, aux_hosp_st_cr = G[g].n_hosp, G[g].n_hosp_star, \
                                                                 G[g].n_hosp_cross, G[g].n_hosp_cr_st
            aux_dead, aux_mild = G[g].n_dead, G[g].n_mild

            aux_omega = 1 - (1 - G[g].omega*dt) * (1 - r0*dt)
##            aux_omega = (G[g].omega + r0)*dt

            ev_inf, ev_inf_acc, ev_hosp, ev_hosp_s, ev_hosp_c, \
                    ev_hosp_cs, ev_mild, ev_dead, flag_beds = md.evol_inf(aux_omega, G[g].xi, G[g].tita, \
                                                                          G[g].tita_star, G[g].tita_cross, G[g].tita_cr_st, \
                                                                          aux_inf0, flag_beds, beds)
            G[g].n_inf += ev_inf
            G[g].n_inf0 += ev_inf
            G[g].n_inf_acc += ev_inf_acc
            G[g].n_hosp += ev_hosp
            G[g].n_hosp_star += ev_hosp_s
            G[g].n_hosp_cross += ev_hosp_c
            G[g].n_hosp_cr_st += ev_hosp_cs
            G[g].n_mild += ev_mild
            G[g].n_dead += ev_dead

            ev_inf, ev_inf_acc, ev_hosp, ev_hosp_s, ev_hosp_c, \
                    ev_hosp_cs, ev_mild, ev_dead, flag_beds = md.evol_inf(aux_omega, G[g].xi, G[g].tita, \
                                                                          G[g].tita_star, G[g].tita_cross, G[g].tita_cr_st, \
                                                                          aux_infA, flag_beds, beds)
            G[g].n_inf += ev_inf
            G[g].n_infA += ev_inf
            G[g].n_inf_acc += ev_inf_acc
            G[g].n_hosp += ev_hosp
            G[g].n_hosp_star += ev_hosp_s
            G[g].n_hosp_cross += ev_hosp_c
            G[g].n_hosp_cr_st += ev_hosp_cs
            G[g].n_mild += ev_mild
            G[g].n_dead += ev_dead

            ev_inf, ev_inf_acc, ev_hosp, ev_hosp_s, ev_hosp_c, \
                    ev_hosp_cs, ev_mild, ev_dead, flag_beds = md.evol_inf(aux_omega, G[g].xi, G[g].tita, \
                                                                          G[g].tita_star, G[g].tita_cross, G[g].tita_cr_st, \
                                                                          aux_infB, flag_beds, beds)
            G[g].n_inf += ev_inf
            G[g].n_infB += ev_inf
            G[g].n_inf_acc += ev_inf_acc
            G[g].n_hosp += ev_hosp
            G[g].n_hosp_star += ev_hosp_s
            G[g].n_hosp_cross += ev_hosp_c
            G[g].n_hosp_cr_st += ev_hosp_cs
            G[g].n_mild += ev_mild
            G[g].n_dead += ev_dead

            #Evolución de asintomáticos
            aux_rec = G[g].n_rec

            ev_mild, ev_asym, ev_inf_acc = md.evol_asym(r0 * dt, G[g].gamma_asym * dt, aux_asym0)
            G[g].n_asym += ev_asym
            G[g].n_asym0 += ev_asym
            G[g].n_mild += ev_mild
            #G[g].n_inf_acc += ev_inf_acc

            ev_mild, ev_asym, ev_inf_acc = md.evol_asym(r0 * dt, G[g].gamma_asym * dt, aux_asymA)
            G[g].n_asym += ev_asym
            G[g].n_asymA += ev_asym
            G[g].n_mild += ev_mild
            #G[g].n_inf_acc += ev_inf_acc

            ev_mild, ev_asym, ev_inf_acc = md.evol_asym(r0 * dt, G[g].gamma_asym * dt, aux_asymB)
            G[g].n_asym += ev_asym
            G[g].n_asymB += ev_asym
            G[g].n_mild += ev_mild
            #G[g].n_inf_acc += ev_inf_acc       

            #Evolución de hospitalizados
            ev_rec, ev_hosp = md.evol_hosp(G[g].gamma * dt, aux_hosp)
            G[g].n_hosp += ev_hosp
            G[g].n_rec += ev_rec

            ev_hosp, ev_hosp_s = md.evol_hosp(G[g].gamma_star * dt, aux_hosp_st)
            G[g].n_hosp_star += ev_hosp_s
            G[g].n_hosp += ev_rec
            flag_beds += ev_hosp_s

            ev_dead, ev_hosp_c = md.evol_hosp(G[g].delta * dt, aux_hosp_cr)
            G[g].n_hosp_cross += ev_hosp_c
            G[g].n_dead += ev_dead

            ev_dead, ev_hosp_cs = md.evol_hosp(G[g].delta_star * dt, aux_hosp_st_cr)
            G[g].n_hosp_cr_st += ev_hosp_cs
            G[g].n_dead += ev_dead
            flag_beds += ev_hosp_cs

            #Evolución de milds
            ev_rec, ev_mild = md.evol_mild(G[g].gamma_mild * dt, aux_mild)
            G[g].n_mild += ev_mild
            G[g].n_rec += ev_rec

            flag_stop += G[g].n_sus + G[g].n_rec + G[g].n_dead
                
        time += 1

    #Fin de la realización
    file.close()
    fileED.close()
    fileICU.close()
    fileSus.close()
    fileAsi.close()
    fileAcu.close()
        
