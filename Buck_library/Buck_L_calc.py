# -*- coding: utf-8 -*-
"""
Created on Wed May  7 10:58:47 2025

@author: daniel.magalhaes_xel
"""

#library for buck indutance calculation
import numpy as np
import Buck_library as bck

def L_calc_1(Vout_max, VDiodo, Dmax, DeltaIL_max, F_max):
    return round(((Vout_max - VDiodo)*(1-Dmax))/(DeltaIL_max*F_max),9)

def AC_DC_ratio_calc (Iout_max, DeltaIL_max):
    return round((DeltaIL_max/Iout_max),2)

def L_calc_2(Vin_min, Vsw, Vout_max, VDiodo, r, F_max, Iout_max):
    return round(((Vin_min-Vsw-Vout_max)*(Vout_max+VDiodo))/(
            (Vin_min-Vsw+VDiodo)*r*F_max*Iout_max),9)

def L_ET_calc (Ton, Vin , Vout, Vsw):
    #infinite number of regulators with different combinations of input and 
    #output voltages but having the same voltseconds are actually the same 
    #regulator from the viewpoint of basic magnetics design. Et is what 
    #really counts. (The only exception to this is the Core Loss term since 
    #this depends directly on the absolute value of the frequency too, not just the Et).
    return round ((Vin - Vsw - Vout) *  Ton,9) #Vusecs

def L_calc_3 (ET_L, I_ratio, Iout):
    return round ( ET_L/(I_ratio * Iout), 9)
 

#final function that will give the range value for the conceptual inductor 
def L_calc (IRatio, Ton, Vin , Vout, Vsw, VDiodo, I_delta_max, F_max, Iout_max):
    Lvalues = np.zeros((1, 3), dtype=float)
    ET_L = bck.L_ET_calc (Ton, Vin, Vout, Vsw)
    Lvalues[0, 0] =  bck.L_calc_1(Vout, VDiodo,Ton, I_delta_max, F_max)
    Lvalues[0, 1] =  bck.L_calc_2(Vin, Vsw, Vout, VDiodo, IRatio, F_max, Iout_max)
    Lvalues[0, 2] = bck.L_calc_3 (ET_L, IRatio, Iout_max)
    return np.max(Lvalues)


def estimation_current_component_AC (ET_L, L):
    return round (ET_L/L, 2)

def estimation_Iratio_calc(Iac_estimated, Iout_max):
    return round (Iac_estimated/Iout_max,2)

def estimation_IL_peak_calc(Iout,Iac_estimated):
    return round (Iout+Iac_estimated/2,2)

def estimation_ILRMS_calc(Iout,Iac_estimated): 
    return round(((Iout**2)+(Iac_estimated**2)/12)**0.5,2)


#-------------------------------------------------------------------

def Estimation_L_copperloss_calc (IL_RMS,L_DCR):
    return round (L_DCR * IL_RMS**2)

def AC_B_Field_calc_1 (L_ET, L_ET100):
    return round (2*(L_ET/L_ET100),2)

def AC_B_Field_calc_2 (L_ET, L_Number_turns, effectiveCoreArea): # using the phisical characteristics
    #effectiveCoreArea is in cm**2
    return round ((100*L_ET/L_Number_turns*effectiveCoreArea),2)

def DC_B_Field_calc_1 (AC_B_Field,Iac_estimated,Iout):
    return round (Iout*AC_B_Field/Iac_estimated,2)

def Peak_B_Field_calc (L_DCB,L_ACB):
    return round (L_DCB+L_ACB/2,2)

def L_energy_core_calc (Ipeak, Lmax):
    return round (0.5*Lmax*Ipeak**2, 9)

def L_Area_product_calc(L_energy, LBm, LJ, LKu) :
    #L_energy - The energy-handling capability of a core
    #LBm - is the flux density, tesla.
    #LJ - is the current density, amps-per-cm2 which controls
    #the copper loss can be seen.
    #LKu - window utilization factor. (which defines the 
    #maximum space that may be used by the copper in the window),
    #LAp - Area product
    LAp = ((2*L_energy)*1e4)/(LBm * LJ * LKu) 
    return round (LAp,9)
    
#def Inductor_core_loss_calc ():