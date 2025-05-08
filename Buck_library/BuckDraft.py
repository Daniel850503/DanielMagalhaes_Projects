"""
Created on Wed May 8 09:33:32 2025

@author: daniel.magalhaes_xel
"""
import numpy as np
import Buck_library as bck

def ILmin_calc(Iout_max, DeltaIL_max):
    return round (Iout_max - ((DeltaIL_max)/2) )

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

def Conceptual_Buck_draft(DCDC_Buck):
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
    return DCDC_Buck