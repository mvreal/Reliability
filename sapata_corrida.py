# -*- coding: utf-8 -*-
"""
Confiabilidade de Sapata Corrida

"""
import numpy as np
from realpy import *

#
def gfunction(x, d):
    c = x[0]
    phi = x[1]
    q0 = x[2]
    gamma = x[3]
    P = x[4]
    B = d[0]
    # Parâmetros de projeto
    Nq = np.exp(np.pi * np.tan(np.radians(phi))) * (np.tan(np.radians(45 + phi/2)))**2
    if phi == 0.00:
        Nc = 2.00 + np.pi
    else:
        Nc = (Nq - 1) * 1 / np.tan(np.radians(phi))
    Ngamma = 2 * (Nq + 1) * np.tan(np.radians(phi))
    # Função de estado limite g(x)=R - S = 0
    r = c * Nc + q0 * Nq + gamma * B / 2 * Ngamma
    q = P / B
    gx = r - q  # Função estado limite
    return gx



xvar = [
    {'varname': 'c', 'vardist': 'normal', 'varmean': 5.00, 'varstd': 1.00},
    {'varname': 'phi', 'vardist': 'normal', 'varmean': 40.00, 'varstd': 4.00},
    {'varname': 'q0', 'vardist': 'normal', 'varmean': 18.00, 'varstd': 2.00},
    {'varname': 'gamma', 'vardist': 'normal', 'varmean': 18.00, 'varstd': 2.00},
    {'varname': 'P', 'vardist': 'normal', 'varmean': 800, 'varstd': 80.00}
    ]

dvar = [{'varname': 'B', 'varvalue': 1.00}]
#
# FORM method
#
sapata_corrida = Reliability(xvar, dvar, gfunction, None, None)
sapata_corrida.form(iHLRF=True, toler=1.e-6)
#