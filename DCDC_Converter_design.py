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
#from Buck_library.Buck_L_calc import Lmax_calc_1 
# Global constants (em UPPER_CASE)


# Configuração de logging (opcional)
"""logging.basicConfig(level=logging.INFO)"""


# Class
class Converter_dcdc :
    def __init__(self, Vin_min, Vin_max, Vout_min, Vout_max, 
                 P_max, F_max, IL_ripple, Vout_ripple, Effmin):
            self.Vin_min = Vin_min
            self.Vin_max = Vin_max
            self.Vout_min = Vout_min
            self.Vout_max = Vout_max
            self.P_max = P_max
            self.F_max = F_max
            self.IL_ripple = IL_ripple # em percentagem
            self.Vout_ripple = Vout_ripple
            self.Effmin = Effmin
            
    def Ganho(self):
        self.Ganho_max = self.Vout_max/self.Vin_min
    
    def Out_current_dc (self):
        self.Iout_dc = round((self.P_max/self.Vout_min),2)
        
    def Out_current_max (self):    
        self.Iout_max = round((self.P_max/self.Vout_min)*
                               (1+((self.IL_ripple/100)/2)), 2)
     
    
    def Delta_current_max (self):
        self.DeltaIL_max = round((self.P_max/self.Vout_min)* 
                                 (self.IL_ripple/100), 2)
     
    def Out_current_min(self):
        self.Iout_min = round (self.Iout_max-self.DeltaIL_max,2)
        
        
    def total_Rdc_Loss_calc( self):
        self.Total_R_loss = round ((self.P_max * ( 1 - self.Effmin / 100 ))/
                                   (self.Iout_max*self.Iout_max))
    
    def total_Vdc_Loss_calc( self ):

        self.Total_V_loss = round ((self.P_max * ( 1 - self.Effmin / 100 ))/
                                   self.Iout_max, 2)
        
    def Loss_Max_calc (self):
        self.total_loss_estimation = round (self.P_max * 
                                    ( 1 - self.Effmin / 100 ),2)
    
    def duty_cycle_estimation_calc (self):
        self.duty_cycle_max = round((self.Vout_max + self.Total_V_loss) / 
                                    self.Vin_min, 2)
        
    def AC_DC_ratio_calc (self):
        self.current_ratio = round((self.DeltaIL_max / self.Iout_dc),2) 
    
    def Periode_calc (self):
        self.T_total = 1/self.F_max
    
    def Switch_time_on_calc (self):
        self.Ton = round (self.T_total * self.duty_cycle_max,9)
    
    def Switch_time_off_calc (self):
        self.Toff = round (self.T_total * (1-self.duty_cycle_max),9)
    
    def Current_ripple_Incrise_rate (self):
        self.Iinc_rate = self.DeltaIL_max / self.Ton
        
    def Current_ripple_decrise_rate (self):
        self.Idec_rate = self.DeltaIL_max / self.Toff

#-------Functions---------


def first_current_converter_simulation ( 
        PWM, Imin,period, increment_on, decrement_off, duty, simulation_step):
    length = PWM.shape[1]
    IL = np.zeros((2, length), dtype=float)
    IL[0] = PWM[0]
    IL[1,0] = Imin
    #step_on = increment_on/(length*duty)
    #step_off = decrement_off/(length*(1-duty))
    step_on = increment_on*simulation_step
    step_off = decrement_off*simulation_step
    for t in range(1, length, 1):
        if (PWM[1,t] == 1):
            IL[1,t] = IL[1,t-1] + step_on
        else:
            IL[1,t] = IL[1,t-1] - step_off
    return IL
    
def first_CAPout_current_converter_simulation ( 
        PWM, Idelta,period, increment_on, decrement_off, duty, simulation_step):
    length = PWM.shape[1]
    IC = np.zeros((2, length), dtype=float)
    IC[0] = PWM[0]
    IC[1,0] = -Idelta/2
    #step_on = increment_on/(length*duty)
    #step_off = decrement_off/(length*(1-duty))
    step_on = increment_on*simulation_step
    step_off = decrement_off*simulation_step
    for t in range(1, length, 1):
        if (PWM[1,t] == 1):
            IC[1,t] = IC[1,t-1] + step_on
        else:
            IC[1,t] = IC[1,t-1] - step_off
    return IC

def first_voltage_inductor_simulation ( 
        PWM, Vin, Vloss, Vout):
    length = PWM.shape[1]
    VL = np.zeros((2, length), dtype=float)
    VL[0] = PWM[0]
    VL[1] = PWM[1]*Vin - (Vloss + Vout)
    # for t in range(0, length, 1):
    #     if (PWM[1,t] == 1):
    #         VL[1,t] = Vin - (Vloss + Vout)
    #     else:
    #         VL[1,t] = - (Vloss + Vout)
    return VL

def first_voltage_capacitor_simulation (IC, Periode, Vout_ripple, Vout_max):

    Up_slope = 0.0
    down_slope = 0.0
    zero_slope = 0.0
    Vout_peak_max = Vout_max + Vout_ripple/2
    Vout_peak_min = Vout_max - Vout_ripple/2
    length = IC.shape[1]
    VC = np.zeros((2, length), dtype=float)
    for i in range (1,length,1):
        #step 1: only for the periode
        if (IC[0,i] <= Periode):
            #step 2: define the ascending slope 
            if (IC[1,i-1] < IC[1,i]):
                Up_slope = Up_slope + 1.0
            #step 2: define the descending slope 
            elif (IC[1,i-1] > IC[1,i]):
                down_slope = down_slope + 1.0
            else: 
                zero_slope = zero_slope + 1.0
        else :
            i = length
        
    Up_slope = Vout_ripple / Up_slope
    down_slope = Vout_ripple / down_slope       
    VC[1,0] = Vout_peak_min
    print("Vout_peak_min = ", VC[1,0])   
    for i in range (1,length,1):
        VC[0,i] = IC[0,i]
        if (IC[1,i-1] < IC[1,i]):
            VC[1,i] =  VC[1,i-1] + Up_slope           
        elif (IC[1,i-1] > IC[1,i]):
            VC[1,i] = VC[1,i-1] - down_slope
        else:
            #não mudar nada
            VC[1,i] = VC[1,i-1]
    return VC


def Duty_cycle_max_buck_calc(Vout_max, VDiodo, Vsw,Vin_min, Effmin):
    return round(((Vout_max+VDiodo)/(Vin_min-Vsw+VDiodo))*(1+(100-Effmin)/100),2)   

def Duty_cycle_max_buck_calc_2(Vout_max, VDiodo, Vsw,Vin_min, Vdc_Loss):
    if VDiodo > Vsw :
        Vsemicondutor = VDiodo
    else:
        Vsemicondutor = Vsw
    return round(((Vout_max+Vdc_Loss+Vsemicondutor)/(Vin_min)),2)  




def ILmin_calc(Iout_max, DeltaIL_max):
    return round (Iout_max - ((DeltaIL_max)/2) )

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
    DCDC_Buck= Converter_dcdc(12.0, 12.0, 1.8, 1.8, 120.0, 500000.0, 30.0, 0.010, 92.0)
    DCDC_Buck.Ganho()
    DCDC_Buck.Delta_current_max()
    DCDC_Buck.Out_current_dc()
    DCDC_Buck.Out_current_max()
    DCDC_Buck.Out_current_min()
    DCDC_Buck.total_Rdc_Loss_calc()
    DCDC_Buck.total_Vdc_Loss_calc()
    DCDC_Buck.duty_cycle_estimation_calc()
    DCDC_Buck.AC_DC_ratio_calc()
    DCDC_Buck.Loss_Max_calc()
    DCDC_Buck.Periode_calc ()
    DCDC_Buck.Switch_time_on_calc()
    DCDC_Buck.Switch_time_off_calc ()
    DCDC_Buck.Current_ripple_Incrise_rate ()   
    DCDC_Buck.Current_ripple_decrise_rate ()
    
    simulation_step = sgn.step_simulation_calc(DCDC_Buck.T_total, 2938)
        
    triangular_wave = sgn.triangular_wave_generation(simulation_step, 
                      DCDC_Buck.T_total, 3*DCDC_Buck.T_total,2)
    
    PWM_signal = sgn.PWM_generation (triangular_wave, DCDC_Buck.duty_cycle_max)
    
    IL = first_current_converter_simulation ( 
            PWM_signal, DCDC_Buck.Iout_min,DCDC_Buck.T_total, 
            DCDC_Buck.Iinc_rate, DCDC_Buck.Idec_rate, DCDC_Buck.duty_cycle_max,
            simulation_step)
    
    IC = first_CAPout_current_converter_simulation( 
            PWM_signal, DCDC_Buck.DeltaIL_max,DCDC_Buck.T_total, 
            DCDC_Buck.Iinc_rate, DCDC_Buck.Idec_rate, DCDC_Buck.duty_cycle_max,
            simulation_step)
    
    VL = first_voltage_inductor_simulation (PWM_signal, DCDC_Buck.Vin_max, 
                                            DCDC_Buck.Total_V_loss, 
                                            DCDC_Buck.Vout_max)
    
    VC = first_voltage_capacitor_simulation (IC, DCDC_Buck.T_total, 
                                             DCDC_Buck.Vout_ripple, 
                                             DCDC_Buck.Vout_max)
    
    
  
    
    ILmin = ILmin_calc(DCDC_Buck.Iout_max, DCDC_Buck.DeltaIL_max)
    
    #Dmax = Duty_cycle_max_buck_calc_2(1.8,VDiodo, Vsw, 12, Vdc_loss)
    
#--------------------------- L calculation -----------------------------------    
    
    Lmax_calc1 =  bck.Lmax_calc_1(1.8, VDiodo, 0.15, 
                             20, 500000)
    
    Iratio = bck.AC_DC_ratio_calc (DCDC_Buck.Iout_max, DCDC_Buck.DeltaIL_max)
    
    Lmax_calc2 = bck.Lmax_calc_2(12, Vsw, 1.8, 
                             VDiodo, Iratio, 500000, 76.6)
    
    ET_L = bck.L_ET_calc (DCDC_Buck.Ton, DCDC_Buck.Vin_min, 
                               DCDC_Buck.Vout_max, Vsw)
    
    Lmax_calc3 = bck.Lmax_calc_3 (ET_L, DCDC_Buck.current_ratio, DCDC_Buck.Iout_max)
    
    #----- Para usar com a bobine escolhida e verificar se cumpre com os requisitos
    Iac_estimated = bck.estimation_current_component_AC(ET_L,Lmax_calc2)   
    
    Iratio_estimated = bck.estimation_Iratio_calc(Iac_estimated,DCDC_Buck.Iout_max)
    
    IL_peak_estimated = bck.estimation_IL_peak_calc(DCDC_Buck.Iout_max,Iac_estimated)
    
    IL_RMS = bck.estimation_ILRMS_calc(DCDC_Buck.Iout_dc,Iac_estimated)
    
    L_copperloss_estimation =  bck.Estimation_L_copperloss_calc (IL_RMS,0.02)
    
    Lenergy = bck.L_energy_core_calc (IL_peak_estimated, Lmax_calc2)
    
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
 
    
  
# Execução direta
if __name__ == "__main__":
    main()