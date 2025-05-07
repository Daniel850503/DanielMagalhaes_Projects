# -*- coding: utf-8 -*-
"""
Created on Wed May  7 13:33:32 2025

@author: daniel.magalhaes_xel
"""

def Cout_max_calc_1(DeltaIL_max, Dmax, Vout_ripple, F_max):
    return round ((DeltaIL_max*0.5 * 0.5*(1-Dmax))/(Vout_ripple * F_max),13)

def Cout_max_calc_2(Dmax, Vout_ripple, Vout_max, Lmax, F_max):
    return round (((1-Dmax)/((Vout_ripple/Vout_max)*8*Lmax*(F_max*F_max))),13)

def Cout_max_calc_3(DeltaIL_max, Ganho_V, Vout_ripple, F_max):
    return round ((DeltaIL_max*(Ganho_V)/(Vout_ripple * F_max)),13)