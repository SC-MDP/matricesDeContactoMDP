import networkx as nx
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
batch = 10
n_days = 300
dt = 1/batch
Ntotal = 0
q0 = 0.0 #leer de file
tass = 0 #leer de file
beds  = 60
r0 = 0.0 #leer de file
aux_r0 = 0.0
eps_a = [0.333 if i < 4 else 0.333 if 4 <= i < 18 else 1 for i in range(21)]

#Lectura parametros epidemicos por rango etario
sim_list = np.loadtxt('ParamEtar/sim_list.txt', comments='#')
sim_num = len(sim_list)
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
c_weights = [0.19, 0.18, 0.30, 0.33] #Fumanelli. Caso pre cuarentena.

for i in range(len(c_list)): #Suma pesada de las matrices
    c_matrixOut = np.add(c_matrixOut, c_weights[i] * c_list[i])

#LOOP LISTA DE SIMULACIONES
print('sim, tass, q0, r0, pico_inf, sat_icu, infacc, deaths, deaths_ev')
if not os.path.exists('res_meta_det'):
    os.makedirs('res_meta_det')
metafile = open('res_meta_det/fasepicos.dat','w')
metafile.write('sim, tass, q0, r0, pico_inf' + '\n')
metafile2 = open('res_meta_det/fasesat.dat','w')
metafile2.write('sim, tass, q0, r0, sat_icu' + '\n')
metafile3 = open('res_meta_det/faseinfac.dat','w')
metafile3.write('sim, tass, q0, r0, ac_inf' + '\n')
metafile3w = open('res_meta_det/faseinfacW.dat','w')
metafile3w.write('sim, tass, q0, r0, ac_infW' + '\n')
metafile3t = open('res_meta_det/faseinfacT.dat','w')
metafile3t.write('sim, tass, q0, r0, ac_infT' + '\n')
metafile3r = open('res_meta_det/faseinfacR.dat','w')
metafile3r.write('sim, tass, q0, r0, ac_infR' + '\n')
metafile4 = open('res_meta_det/fasedeath.dat','w')
metafile4.write('sim, tass, q0, r0, death' + '\n')
metafile5 = open('res_meta_det/fasedeathev.dat','w')
metafile5.write('sim, tass, q0, r0, death_ev' + '\n')
metafile6 = open('res_meta_det/fasenr0.dat','w')
metafile6.write('sim, tass, q0, r0, n_tests' + '\n')
metafile7 = open('res_meta_det/n_beds.dat','w')
metafile7.write('sim, tass, q0, r0, n_beds' + '\n')

for sim in range(sim_num):
   
    #Lista de redes. Cada red es un grupo etario
    G = []
    Ntotal = 0
    Inf_init = 0.0

    time = 0
    pico_inf = 0
    dead_ev = 0
    sat_icu = 0
    n_tests = 0.0
    n_tests_pos = 0.0
    n_beds = 0.0
    n_beds_max = 0.0
    n_beds_old = 0.0
    tass, q0, r0 = sim_list[sim]
    aux_r0 = r0

    flag_beds_det = 0
    flag_beds_count = 0

    for i in range(grupos_etar):
        G.append(Grupo())
        G[i].N = pob[i] * 1.060529347048 #Actualizacion a 2020
        
        #Subdivisión del grupo. "0" son los esenciales, A y B son los grupos que se alternan
        G[i].NumA = G[i].N * (1 - q0) / 2
        G[i].NumB = G[i].N * (1 - q0) / 2
        G[i].Num0 = G[i].N - G[i].NumA - G[i].NumB
     
        #Parametros epidémicos del grupo
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
        G[i].S = S_input[i]
        G[i].E = E_input[i]
        G[i].I = I_input[i]
        G[i].A = A_input[i]
        G[i].H =  HnL_input[i]
        G[i].Hstar = HsL_input[i]
        G[i].Hcross = HnD_input[i]
        G[i].Hcrst = HsD_input[i]
        G[i].M = M_input[i]
        G[i].R = R_input[i]
        G[i].D = D_input[i]
        
        G[i].Iacc = Iacc_input/grupos_etar #No importa distribucion porque sacamos la suma
        G[i].Iacc_w = Iacc_input/grupos_etar
        G[i].Iacc_t = 0
        G[i].Iacc_real = 0 #R_input[i] + D_input[i]
        G[i].hICU = G[i].Hstar + G[i].Hcrst

        G[i].Hstar_binf = G[i].Hstar
        G[i].Hcrst_binf = G[i].Hcrst
        G[i].hICU_binf = G[i].Hstar_binf + G[i].Hcrst_binf
        G[i].n_tests = 0
        G[i].n_tests_pos = 0

        #Variable auxiliar
        Ntotal += G[i].N
        Inf_init += G[i].I
        flag_beds_det += G[i].Hstar + G[i].Hcrst
        n_beds += G[i].Hstar_binf + G[i].Hcrst_binf    

    #Estados extra por estrategia
    for g in range(grupos_etar):
        G[g].SA = G[g].S * (1 - q0) / 2
        G[g].SB = G[g].S * (1 - q0) / 2
        G[g].S0 = G[g].S - G[g].SA - G[g].SB

        G[g].EA = G[g].E * (1 - q0) / 2
        G[g].EB = G[g].E * (1 - q0) / 2
        G[g].E0 = G[g].E - G[g].EA - G[g].EB

        G[g].IA = G[g].I * (1 - q0) / 2
        G[g].IB = G[g].I * (1 - q0) / 2
        G[g].I0 = G[g].I - G[g].IA - G[g].IB

        G[g].AA = G[g].A * (1 - q0) / 2
        G[g].AB = G[g].A * (1 - q0) / 2
        G[g].A0 = G[g].A - G[g].AA - G[g].AB

        G[g].f0 = G[g].fA = G[g].fB = 0

    #Output
    ICUtot = 0
    ICUtot_binf = 0
    Itot = 0  
    Dtot = 0      
    Sustot = 0
    Asitot = 0
    Acutot = 0
    Acutot_w = 0
    Acutot_t = 0
    Acutot_real = 0
    Rectot = 0

##    #Creo archivo.
    ffnn2 = '_q0'+str(format(q0, '03g'))+'_tao'+str(tass)+'_r'+str(format(aux_r0, '03g'))
    if not os.path.exists('results_det%s' %ffnn2):
        os.makedirs('results_det%s' %ffnn2)
##        
    file2=open('results_det%s/DetInf.dat' %ffnn2,'w')
    fileICU=open('results_det%s/DetICU.dat' %ffnn2,'w')
    fileICUbinf=open('results_det%s/DetICUbinf.dat' %ffnn2,'w')
    fileSus=open('results_det%s/DetSus.dat' %ffnn2,'w')
    fileD=open('results_det%s/DetDead.dat' %ffnn2,'w')
    fileAsi=open('results_det%s/DetAsi.dat' %ffnn2,'w')
    fileAcu=open('results_det%s/DetAcu.dat' %ffnn2,'w')
    fileAcuW=open('results_det%s/DetAcuW.dat' %ffnn2,'w')
    fileAcuT=open('results_det%s/DetAcuT.dat' %ffnn2,'w')
    fileAcuR=open('results_det%s/DetAcuR.dat' %ffnn2,'w')
    fileRec=open('results_det%s/DetRec.dat' %ffnn2,'w')

    #DINAMICA
    while time < n_days*batch:

        #Cantidades acumuladas       
        Aux=0
        for g in range(grupos_etar):
            Aux = Aux + G[g].hICU
        ICUtot = Aux

        Aux=0
        for g in range(grupos_etar):
            Aux = Aux + G[g].hICU_binf
        ICUtot_binf = Aux

        Aux=0
        for g in range(grupos_etar):
            Aux = Aux + G[g].n_tests
        n_tests = Aux

        Aux=0
        for g in range(grupos_etar):
            Aux = Aux + G[g].n_tests_pos
        n_tests_pos = Aux
        
        Aux=0
        for g in range(grupos_etar):
            Aux = Aux + G[g].I
        Itot = Aux
        
        if Itot > pico_inf:
            pico_inf = Itot #Registro maximo de inf
        
        Aux=0
        for g in range(grupos_etar):
            Aux = Aux + G[g].D
        Dtot = Aux

        Aux=0
        for g in range(grupos_etar):
            Aux = Aux + G[g].S
        Sustot = Aux              

        Aux=0
        for g in range(grupos_etar):
            Aux = Aux + G[g].A
        Asitot = Aux

        Aux=0
        for g in range(grupos_etar):
            Aux = Aux + G[g].Iacc
        Acutot = Aux

        Aux=0
        for g in range(grupos_etar):
            Aux = Aux + G[g].Iacc_w
        Acutot_w = Aux

        Aux=0
        for g in range(grupos_etar):
            Aux = Aux + G[g].Iacc_t
        Acutot_t = Aux

        Aux=0
        for g in range(grupos_etar):
            Aux = Aux + G[g].Iacc_real
        Acutot_real = Aux

        Aux=0
        for g in range(grupos_etar):
            Aux = Aux + G[g].R
        Rectot = Aux
        
        #Salvo en los archivos
        if time%batch==0:
            file2.write(str(time*dt) + "\t" + str(Itot) + "\n")
            fileICU.write(str(time*dt) + "\t" + str(ICUtot) + "\n")
            fileICUbinf.write(str(time*dt) + "\t" + str(ICUtot_binf) + "\n")
            fileD.write(str(time*dt) + "\t" + str(Dtot) + "\n")
            fileSus.write(str(time*dt) + "\t" + str(Sustot) + "\n")
            fileAsi.write(str(time*dt) + "\t" + str(Asitot) + "\n")
            fileAcu.write(str(time*dt) + "\t" + str(Acutot) + "\n")
            fileAcuW.write(str(time*dt) + "\t" + str(Acutot_w) + "\n")
            fileAcuT.write(str(time*dt) + "\t" + str(Acutot_t) + "\n")
            fileAcuR.write(str(time*dt) + "\t" + str(Acutot_real) + "\n")
            fileRec.write(str(time*dt) + "\t" + str(Rectot) + "\n")

        #Estrategia de tests
        if Itot <= Inf_init: #corto estrategia de tests
##        if Itot <= 1:
            r0 = 0.0
        else:
            r0 = aux_r0
        
        #Calculo de las fuerzas de infección
        if math.sin(math.pi*(time*dt/tass))>=0:
            for g in range(grupos_etar):
                G[g].fout = G[g].fin = 0
                    
                for i in range(grupos_etar):
                    G[g].fout += c_matrixOut[g][i]*(G[g].beta*(G[i].I0+G[i].IA)+G[g].beta_a*(G[i].A0+G[i].AA)) / G[i].N #deterministas
                    
                G[g].f0, G[g].fA, G[g].fB = G[g].fout, G[g].fout, G[g].fin 
                
        else:
            for g in range(grupos_etar):
                G[g].fout = G[g].fin = 0
                    
                for i in range(grupos_etar):
                    G[g].fout += c_matrixOut[g][i]*(G[g].beta*(G[i].I0+G[i].IB)+G[g].beta_a*(G[i].A0+G[i].AB)) / G[i].N #deterministas
                    
                G[g].f0, G[g].fB, G[g].fA = G[g].fout, G[g].fout, G[g].fin 
          
        #Evolución de la enfermedad    
        for g in range(grupos_etar):

            #Deterministas
            G[i].n_tests = G[i].n_tests + dt * r0 * (G[g].S + G[g].E + G[g].I + G[g].A)
            G[i].n_tests_pos = G[i].n_tests_pos + dt * r0 * (G[g].I + G[g].A)
            
            aux_S0 = G[g].S0 - dt * G[g].f0 * G[g].S0
            aux_SA = G[g].SA - dt * G[g].fA * G[g].SA
            aux_SB = G[g].SB - dt * G[g].fB * G[g].SB
            G[g].S = aux_S0 + aux_SA + aux_SB
            
            aux_E0 = G[g].E0 + dt * (G[g].f0 * G[g].S0 - G[g].alpha * G[g].E0)
            aux_EA = G[g].EA + dt * (G[g].fA * G[g].SA - G[g].alpha * G[g].EA)
            aux_EB = G[g].EB + dt * (G[g].fB * G[g].SB - G[g].alpha * G[g].EB)    
            G[g].E = aux_E0 + aux_EA + aux_EB

            aux_Iacc = G[g].Iacc + dt * (G[g].omega + r0) * (G[g].I0 + G[g].IA + G[g].IB)
            aux_Iacc_w = G[g].Iacc_w + dt * G[g].omega * (G[g].I0 + G[g].IA + G[g].IB)
            aux_Iacc_t = G[g].Iacc_t + dt * r0 * (G[g].I0 + G[g].IA + G[g].IB) + dt * r0 * G[g].A
            
            G[g].I  = G[g].I0 + G[g].IA + G[g].IB   #estos son los infectados viejos por lo que no hay que poner "aux"
            aux_I0 = G[g].I0  + dt * (G[g].epsilon * G[g].alpha * G[g].E0 - (G[g].omega + r0) * G[g].I0)
            aux_IA = G[g].IA  + dt * (G[g].epsilon * G[g].alpha * G[g].EA - (G[g].omega + r0) * G[g].IA)       
            aux_IB = G[g].IB  + dt * (G[g].epsilon * G[g].alpha * G[g].EB - (G[g].omega + r0) * G[g].IB)   
            aux_I = aux_I0 + aux_IA + aux_IB #No actualizo I aun para usarlo a tiempo t en otras ecs
            
            aux_A0 = G[g].A0 + dt * ((1 - G[g].epsilon) * G[g].alpha * G[g].E0 - G[g].gamma_asym * G[g].A0 - r0 * G[g].A0)      
            aux_AA = G[g].AA + dt * ((1 - G[g].epsilon) * G[g].alpha * G[g].EA - G[g].gamma_asym * G[g].AA - r0 * G[g].AA)
            aux_AB = G[g].AB + dt * ((1 - G[g].epsilon) * G[g].alpha * G[g].EB - G[g].gamma_asym * G[g].AB - r0 * G[g].AB)
            aux_A = aux_A0 + aux_AA + aux_AB

            aux_old_beds = G[g].Hstar + G[g].Hcrst
            aux_H = G[g].H + dt * (G[g].tita * G[g].xi * (G[g].omega + r0) * G[g].I + G[g].gamma_star * G[g].Hstar- G[g].gamma * G[g].H)
            aux_Hstar = G[g].Hstar + dt * (G[g].tita_star * G[g].xi * (G[g].omega + r0) * G[g].I - G[g].gamma_star * G[g].Hstar)
            aux_Hcross = G[g].Hcross + dt * (G[g].tita_cross * G[g].xi * (G[g].omega + r0) * G[g].I - G[g].delta * G[g].Hcross)
            aux_Hcrst = G[g].Hcrst + dt * (G[g].tita_cr_st * G[g].xi * (G[g].omega + r0) * G[g].I - G[g].delta_star * G[g].Hcrst)
            aux_M = G[g].M + dt * ((1 - G[g].xi) * (G[g].omega + r0) * G[g].I - G[g].gamma_mild * G[g].M + r0 * G[g].A)
            aux_R = G[g].R + dt * (G[g].gamma * G[g].H + G[g].gamma_mild * G[g].M + G[g].gamma_asym * G[g].A)
            aux_D = G[g].D + dt * (G[g].delta * G[g].Hcross + G[g].delta_star * G[g].Hcrst)

            aux_RD = aux_R + aux_D

            flag_beds_det += aux_Hstar + aux_Hcrst - aux_old_beds

            n_beds_old = G[g].Hstar_binf + G[g].Hcrst_binf
            aux_Hstar_binf = G[g].Hstar_binf + dt * (G[g].tita_star * G[g].xi * (G[g].omega + r0) * G[g].I - G[g].gamma_star * G[g].Hstar_binf)
            aux_Hcrst_binf = G[g].Hcrst_binf + dt * (G[g].tita_cr_st * G[g].xi * (G[g].omega + r0) * G[g].I - G[g].delta_star * G[g].Hcrst_binf)
            
            n_beds += aux_Hstar_binf + aux_Hcrst_binf - n_beds_old
            if n_beds > n_beds_max:
                n_beds_max = n_beds
                
            if flag_beds_det >= beds:

                sup = flag_beds_det - beds
                dead_ev += sup
                
                delHstar = G[g].tita_star * G[g].xi * (G[g].omega + r0) * G[g].I - G[g].gamma_star * G[g].Hstar
                delHcrst = G[g].tita_cr_st * G[g].xi * (G[g].omega + r0) * G[g].I - G[g].delta_star * G[g].Hcrst
                
                aux_Hstar = aux_Hstar - sup * delHstar / (delHstar + delHcrst)
                aux_Hcrst = aux_Hcrst - sup * delHcrst / (delHstar + delHcrst)
                aux_D = aux_D + sup  
                
                flag_beds_det = beds
                flag_beds_count = 1
                             
            G[g].S0 = aux_S0
            G[g].SA = aux_SA 
            G[g].SB = aux_SB
            G[g].E0 = aux_E0
            G[g].EA = aux_EA
            G[g].EB = aux_EB
            G[g].I0 = aux_I0
            G[g].IA = aux_IA
            G[g].IB = aux_IB
            G[g].A0 = aux_A0
            G[g].AA = aux_AA
            G[g].AB = aux_AB           
            G[g].M = aux_M
            G[g].R = aux_R
            #New block
            G[g].Iacc = aux_Iacc
            G[g].Iacc_w = aux_Iacc_w
            G[g].Iacc_t = aux_Iacc_t
            G[g].Iacc_real = aux_RD
            #----------------------------
            G[g].hICU = aux_Hstar + aux_Hcrst
            G[g].hICU_binf = aux_Hstar_binf + aux_Hcrst_binf
            G[g].I = aux_I
            G[g].H = aux_H
            G[g].Hstar = aux_Hstar
            G[g].Hcross = aux_Hcross
            G[g].Hcrst = aux_Hcrst
            G[g].Hstar_binf = aux_Hstar_binf
            G[g].Hcrst_binf = aux_Hcrst_binf
            G[g].A = aux_A
            G[g].D = aux_D

        time += 1
        
        if flag_beds_count == 1:
            sat_icu += 1 #cuanto tiempo estan saturadas las camas
            flag_beds_count = 0
            

    #Fin de la realización
    file2.close()
    fileICU.close()
    fileICUbinf.close()
    fileSus.close()
    fileD.close()
    fileAsi.close()
    fileAcu.close()
    fileAcuW.close()
    fileAcuT.close()
    fileAcuR.close()
    fileRec.close()

    print(sim, tass, q0, aux_r0, pico_inf, sat_icu/batch, Acutot, Dtot, dead_ev)
    metafile.write(str(sim) + '\t' + str(tass) + '\t' + str(q0) + '\t' + str(aux_r0) + '\t' + str(pico_inf) + '\n')
    metafile2.write(str(sim) + '\t' + str(tass) + '\t' + str(q0) + '\t' + str(aux_r0) + '\t' + str(sat_icu/batch) + '\n')
    metafile3.write(str(sim) + '\t' + str(tass) + '\t' + str(q0) + '\t' + str(aux_r0) + '\t' + str(Acutot) + '\n')
    metafile3w.write(str(sim) + '\t' + str(tass) + '\t' + str(q0) + '\t' + str(aux_r0) + '\t' + str(Acutot_w) + '\n')
    metafile3t.write(str(sim) + '\t' + str(tass) + '\t' + str(q0) + '\t' + str(aux_r0) + '\t' + str(Acutot_t) + '\n')
    metafile3r.write(str(sim) + '\t' + str(tass) + '\t' + str(q0) + '\t' + str(aux_r0) + '\t' + str(Acutot_real) + '\n')
    metafile4.write(str(sim) + '\t' + str(tass) + '\t' + str(q0) + '\t' + str(aux_r0) + '\t' + str(Dtot) + '\n')
    metafile5.write(str(sim) + '\t' + str(tass) + '\t' + str(q0) + '\t' + str(aux_r0) + '\t' + str(dead_ev) + '\n')
    metafile6.write(str(sim) + '\t' + str(tass) + '\t' + str(q0) + '\t' + str(aux_r0) + '\t' + str(n_tests) + '\t' + str(n_tests_pos) + '\n')
    metafile7.write(str(sim) + '\t' + str(tass) + '\t' + str(q0) + '\t' + str(aux_r0) + '\t' + str(n_beds_max) + '\n')

metafile.close()
metafile2.close()
metafile3.close()
metafile3w.close()
metafile3t.close()
metafile3r.close()
metafile4.close()
metafile5.close()
metafile6.close()
metafile7.close()
