# -*- coding: utf-8 -*-
"""
Created on Wed May  7 14:04:04 2025

@author: daniel.magalhaes_xel
"""

import numpy as np

def step_simulation_calc(Periodo, total_samples_number):
    return Periodo/total_samples_number

def triangular_wave_1(step, Periodo, total_time):
    Total_number_samples = int (round(total_time / step , 0))
    sample_triangular = int (round(Periodo / step , 0))
    sample_peak = sample_triangular / 2
    triangular_increment = 1 / sample_peak
    matriz = np.zeros((2, Total_number_samples), dtype=float)
    
    #matriz[0, t] -> tempo 
    #matriz[1, t] -> wave
    Incrise = 1;
    for t in range(1, Total_number_samples, 1):
        matriz[0, t] = round ( matriz[0, t-1] + step, 12)
        if (matriz[1, t-1]< (sample_peak/sample_peak)) and Incrise == 1:
            matriz[1, t] = matriz[1, t-1] + triangular_increment
        elif matriz[1, t-1] <= 0:
            Incrise = 1
            matriz[1, t] = matriz[1, t-1] + triangular_increment
        else:
           Incrise = 0
           matriz[1, t] = matriz[1, t-1] - triangular_increment
    return matriz

def triangular_wave_2(step, Periodo, total_time):
    Total_number_samples = int (round(total_time / step , 0))
    sample_triangular = int (round(Periodo / step , 0))
    sample_peak = sample_triangular
    triangular_increment = 1 / sample_peak
    matriz = np.zeros((2, Total_number_samples), dtype=float)
    
    #matriz[0, t] -> tempo 
    #matriz[1, t] -> wave
    Incrise = 1;
    for t in range(1, Total_number_samples, 1):
        matriz[0, t] = round ( matriz[0, t-1] + step, 12)
        if (matriz[1, t-1]< (sample_peak/sample_peak)) and Incrise == 1:
            matriz[1, t] = matriz[1, t-1] + triangular_increment
        elif matriz[1, t-1] <= 0:
            Incrise = 1
            matriz[1, t] = matriz[1, t-1] + triangular_increment
        else:
           Incrise = 0
           matriz[1, t] = 0
    return matriz

def triangular_wave_generation(step, Periodo, total_time, type_signal):
    if (type_signal == 1):
        matriz = triangular_wave_1(step, Periodo, total_time)
    elif (type_signal == 2):
        matriz = triangular_wave_2(step, Periodo, total_time)
    else:
        print("ERROR IN SIGNAL TYPE DEFINITION. ONLY 1 OR 2")
        Total_number_samples = int (round(total_time / step , 0))
        matriz = np.zeros((2, Total_number_samples), dtype=float)
    return matriz

def PWM_generation (triangular_wave, duty_cycle):
    length = triangular_wave.shape[1]
    PWM = np.zeros((2, length), dtype=float)
    PWM[0] = triangular_wave[0]
    for t in range(0, length, 1):
        if (triangular_wave[1,t] <= duty_cycle):
            PWM[1,t] = 1
        else:
            PWM[1,t] = 0
    return PWM