# -*- coding: utf-8 -*-
"""
Created on Wed May  7 11:21:53 2025

@author: daniel.magalhaes_xel
"""

from .Buck_L_calc import (  
    Lmax_calc_1,
    AC_DC_ratio_calc,
    Lmax_calc_2,
    L_ET_calc,
    Lmax_calc_3,
    estimation_current_component_AC,
    estimation_Iratio_calc,
    estimation_IL_peak_calc,
    estimation_ILRMS_calc,
    Estimation_L_copperloss_calc,
    AC_B_Field_calc_1,
    AC_B_Field_calc_2,
    DC_B_Field_calc_1,
    Peak_B_Field_calc,
    L_energy_core_calc,
    L_Area_product_calc )

from .Buck_C_calc import (Cout_max_calc_1,
    Cout_max_calc_2,
    Cout_max_calc_3)

__all__ = [
    "Lmax_calc_1",
    "AC_DC_ratio_calc",
    "Lmax_calc_2",
    "L_ET_calc",
    "Lmax_calc_3",
    "estimation_current_component_AC",
    "estimation_Iratio_calc",
    "estimation_IL_peak_calc",
    "estimation_ILRMS_calc",
    "Estimation_L_copperloss_calc",
    "AC_B_Field_calc_1",
    "AC_B_Field_calc_2",
    "DC_B_Field_calc_1",
    "Peak_B_Field_calc",
    "L_energy_core_calc",
    "L_Area_product_calc",
    "Cout_max_calc_1",
    "Cout_max_calc_2",
    "Cout_max_calc_3"
]