# -*- coding: utf-8 -*-
"""
Created on Wed May  7 14:28:15 2025

@author: daniel.magalhaes_xel
"""
import matplotlib.pyplot as plt

def plot_func (x, y , titulo_grafico, x_legenda, y_legenga):
    plt.plot(x, y)
    plt.title(titulo_grafico)
    plt.xlabel(x_legenda)
    plt.ylabel(y_legenga)
    plt.grid(True)
    plt.show()
    return 1

def plot_matrix_2xN_func(matriz, Plot_title, x_title, y_title):
    x = matriz [0]
    y = matriz [1]
    plot_func (x, y , Plot_title, x_title, y_title)
