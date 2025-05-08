# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 14:04:59 2025

@author: daniel.magalhaes_xel
"""

# Imports
import os
import math 
import statistics
import matplotlib.pyplot as plt
import numpy as np
import Buck_library as bck
import Signals_library as sgn
import Plot_library as plt

# Global constants (em UPPER_CASE)


# Configuração de logging (opcional)

"""logging.basicConfig(level=logging.INFO)"""

def main():
    print("Dimencionamento do conversor")
    #VARIAVEIS LOCAIS
    VDiodo=0.5 # apenas para teste
    Vsw=1.5 # apenas para teste
    tempo = []
    IL_tempo = []
    Transistor_control = []
    periodo = 0.0
    passo = 0.0
    Ton = 0.0
    Passo_corrente_L_Ton  = 0.0
    Passo_corrente_L_Toff = 0.0
    N_amostras = 110
    Lmax_calc = 0.0
    Idelta_max=0.0
    ET_L=0;
    #CORPO DO PROGRAMA
    DCDC_Buck = bck.Converter_dcdc(12.0, 12.0, 1.8, 1.8, 120.0, 500000.0, 30.0, 0.010, 92.0)
    DCDC_Buck = bck.Conceptual_Buck_draft(DCDC_Buck)

    simulation_step = sgn.step_simulation_calc(DCDC_Buck.T_total, 2938)
        
    triangular_wave = sgn.triangular_wave_generation(simulation_step, 
                      DCDC_Buck.T_total, 3*DCDC_Buck.T_total,2)
    
    PWM_signal = sgn.PWM_generation (triangular_wave, DCDC_Buck.duty_cycle_max)
    
    IL = bck.first_current_converter_simulation ( 
            PWM_signal, DCDC_Buck.Iout_min,DCDC_Buck.T_total, 
            DCDC_Buck.Iinc_rate, DCDC_Buck.Idec_rate, DCDC_Buck.duty_cycle_max,
            simulation_step)
    
    IC = bck.first_CAPout_current_converter_simulation( 
            PWM_signal, DCDC_Buck.DeltaIL_max,DCDC_Buck.T_total, 
            DCDC_Buck.Iinc_rate, DCDC_Buck.Idec_rate, DCDC_Buck.duty_cycle_max,
            simulation_step)
    
    VL = bck.first_voltage_inductor_simulation (PWM_signal, DCDC_Buck.Vin_max, 
                                            DCDC_Buck.Total_V_loss, 
                                            DCDC_Buck.Vout_max)
    
    VC = bck.first_voltage_capacitor_simulation (IC, DCDC_Buck.T_total, 
                                             DCDC_Buck.Vout_ripple, 
                                             DCDC_Buck.Vout_max)
    
    
  
    
    ILmin = bck.ILmin_calc(DCDC_Buck.Iout_max, DCDC_Buck.DeltaIL_max)
    
 

    #Dmax = Duty_cycle_max_buck_calc_2(1.8,VDiodo, Vsw, 12, Vdc_loss)
    
#--------------------------- L calculation -----------------------------------        
    L =bck.L_calc (DCDC_Buck.current_ratio, DCDC_Buck.Ton, DCDC_Buck.Vin_min,
               DCDC_Buck.Vout_max, Vsw, VDiodo, DCDC_Buck.DeltaIL_max, 
               DCDC_Buck.F_max, DCDC_Buck.Iout_max)
    
    print("Lmin: ", L, "H")

    """
    #----- Para usar com a bobine escolhida e verificar se cumpre com os requisitos
    Iac_estimated = bck.estimation_current_component_AC(ET_L,Lmax_calc2)   
    
    #Iratio_estimated = bck.estimation_Iratio_calc(Iac_estimated,DCDC_Buck.Iout_max)
    
    #IL_peak_estimated = bck.estimation_IL_peak_calc(DCDC_Buck.Iout_max,Iac_estimated)
    
    #IL_RMS = bck.estimation_ILRMS_calc(DCDC_Buck.Iout_dc,Iac_estimated)
    
    #L_copperloss_estimation =  bck.Estimation_L_copperloss_calc (IL_RMS,0.02)
    
    #Lenergy = bck.L_energy_core_calc (IL_peak_estimated, Lmax_calc2)
    
    Cout_min_1 = bck.Cout_max_calc_1(20, 0.165, 
                                 0.01, 600000)
   
    Cout_min_2 = bck.Cout_max_calc_2(0.165,0.01, 
                                 1.8, 150E-9,600000)
    
    Cout_min_3 = bck.Cout_max_calc_3 ( 20, 0.15 ,
                                  0.01, 600000 )
    
    
    print (" ------------------- ESTIMAÇÃO INICIAL ------------------------ ")
    print ("-----------------------------------------------------------------")
    print ("Tensão de entrada minima: ", DCDC_Buck.Vin_min, "V")
    print ("Tensão de entrada minima: ", DCDC_Buck.Vout_max, "V")
    print ("Ganho maximo:", DCDC_Buck.Ganho_max)
    print ("Perdas totais maximas esperadas :", 
          DCDC_Buck.total_loss_estimation, "V")
    print ("Tensão equivalente às Perdas totais maximas esperadas :", 
          DCDC_Buck.Total_V_loss, "V")
    print ("Corrente DC de saída", DCDC_Buck.Iout_dc, "A")
    print ("Corrente de pico de saída", DCDC_Buck.Iout_max, "A")
    print ("Corrente de pico de saída", DCDC_Buck.Iout_min, "A")
    print ("variação maxima da Corrente de saída", DCDC_Buck.DeltaIL_max, "A")
    print ("Duty cycle max estimado: ", DCDC_Buck.duty_cycle_max)
    print ("Tempo do transistor em condução: ", DCDC_Buck.Ton, "s")
    print ("Tempo do transistor em condução: ", DCDC_Buck.Toff, "s")
    print ("Relação entre componente AC e DC da corrente:", Iratio)
    print ("Declive da rampa de subida da currente: ", 
           DCDC_Buck.Iinc_rate*1e-6, "A/us")   
    print ("Declive da rampa de descida da currente: ", 
           DCDC_Buck.Idec_rate*1e-6, "A/us")
    print (" ---------------- first simulation Data----------------")
    print ("Total de amostras da simulação: ", triangular_wave[0,100] )
    print ("-------------------- CALCULO DA BOBINE --------------------------")
    print ("Bobine minima necessaria calculo 1: ",  Lmax_calc1, "H")
    print ("Bobine minima necessaria calculo 2: ",  Lmax_calc2, "H")
    print ("Energia da Bobine - voltsecs: ",  ET_L, "Voltsec")
    print ("Bobine minima necessaria calculo 3: ",  Lmax_calc3, "H")
    print ("Eneria da Bobine: ",  Lenergy, "J")
    print ("------------------- Dinamic estimation --------------------------")
    print ("Componente AC da corrente: ",Iac_estimated,"A")
    print ("Estimação da relação entre componente AC e DC da corrente:"
           , Iratio_estimated)
    print ("Estimação da currente eficaz da bobine: ", IL_RMS, "A")
    print ("-------------------- CALCULO DO CONDENSADOR DE SAÍDA ------------")
    print ("condensador de saída minimo - calculo1: ", Cout_min_1, "F")
    print ("condensador de saída minimo - calculo2: ", Cout_min_2, "F")
    print ("condensador de saída minimo - calculo3: ", Cout_min_3, "F")
    print ("-----------------------------------------------------------------")
    
    imagem1 = plt.plot_matrix_2xN_func(triangular_wave, "triangular", "tempo", "sinal")
    imagem2 = plt.plot_matrix_2xN_func(PWM_signal, "PWM", "tempo", "sinal")
    imagem3 = plt.plot_matrix_2xN_func(IL, "Inductor current", "time", "Current")
    imagem4 = plt.plot_matrix_2xN_func(IC, "Capacitor current", "time", "Current")
    imagem5 = plt.plot_matrix_2xN_func(VL, "Inductor voltage", "time", "Voltage")
    imagem6 = plt.plot_matrix_2xN_func(VC, "Capacitor voltage", "time", "Voltage")


    # Dmax = 0.1745
    # Lmax_calc = Lmax_calc2
    # periodo = 1.0/(500000.0)
    # passo= periodo / (N_amostras-1)
    
    # Ton = Dmax*periodo
    # Idelta_max = ((1/Lmax_calc)*
    #              (DCDC_Buck.Vin_min-Vsw-Vdc_loss-DCDC_Buck.Vout_max))*Ton
    
    # print ("Idelta_max: ", Idelta_max)
    # Passo_corrente_L_Ton = (Idelta_max / (periodo*Dmax)) * passo
    # Passo_corrente_L_Toff = (Idelta_max / (periodo*(1-Dmax))) * passo
                 
    # print("current minima da bobine em regime permanente: ", ILmin)
    # print( "passo da corrente ON: ", Passo_corrente_L_Ton)
    # print( "passo da corrente OFF: ", Passo_corrente_L_Toff)
    # print("Ton: ", Ton)
    
    # tempo.append(0) 
    # Transistor_control.append(0) 
    # IL_tempo.append(ILmin)
         
    # for t in range(1, N_amostras, 1):
    # #     print(" posição do array: ", t-1)
    #     tempo.append(tempo[t-1]+passo)
    # #     #print("tempo: ",tempo[t-1])
    #     if tempo[t-1] < ( Dmax * periodo ):
    #          #print ("está em ton: ",tempo[t-1],"s" )
    #           Transistor_control.append(1) 
    #           IL_tempo.append(IL_tempo[t-1]+Passo_corrente_L_Ton)
    #     else:
    #          Transistor_control.append(0)
    #          IL_tempo.append(IL_tempo[t-1]-Passo_corrente_L_Toff)
    
    # tempo.append(periodo)
    # Transistor_control.append(1)    
    # IL_tempo.append(ILmin)
    
    # print("pico de corrente estimado:", max(IL_tempo), "A" )
    # print("média de corrente estimado:",  statistics.mean(IL_tempo), "A" )
    # plot1 = plot_func (tempo, IL_tempo, "Currente na bobime", "tempo [s]", "currente [A]")
    # plot2 = plot_func (tempo, Transistor_control, "controlo do transistor", 
    #                    "tempo [s]", "estado do transistor")
 
    
""" 
# Execução direta
if __name__ == "__main__":
    main()