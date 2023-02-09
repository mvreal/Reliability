# -*- coding: utf-8 -*-
"""
Reliability of reinforced concrete beams in fire situation

"""
import numpy as np
from realpy import *

#
def gfunction(x, d):
    P = x[0]
    H = x[1]
    E = x[2]
    L = x[3]
    I = x[4]
    b = x[5]
   
    # Parâmetro limite de deslocamento superior
    delta_lim = d[0]  # Deslocamento limite  no topo da coluna
    
    
    # Função de estado limite g(x)=delta_lim - delta = 0
    psi = L * np.sqrt(P/(E*I))
    t = H / (E*I*(P/(E*I))**1.5)
    a0 = np.tan(psi) - psi
    a1 = H*E*I*(np.tan(psi))**2
    a2 = 4*H*E*I*np.tan(psi)
    num = 1-4*(a1/b**2)
    if num<0.00: 
        a3 = 0.00
    else: 
        a3 = num
    a4 = 1 + np.sqrt(a3)
  
    delta = t*(a0 + b**2*a4**2/a2)

    

    gx = delta_lim - delta  # Função estado limite
    return gx





xvar = [
    {'varname': 'P', 'vardist': 'normal', 'varmean': 1000., 'varcov': 0.2},
    {'varname': 'H', 'vardist': 'normal', 'varmean': 10., 'varcov': 0.2},
    {'varname': 'E', 'vardist': 'normal', 'varmean': 20.e6, 'varcov': 0.1},
    {'varname': 'L', 'vardist': 'normal', 'varmean': 10., 'varcov': 0.05},
    {'varname': 'I', 'vardist': 'normal', 'varmean': 1.0417e-2, 'varcov': 0.05},
    {'varname': 'b', 'vardist': 'normal', 'varmean': 0.4167e-6, 'varcov': 0.1}
    
]

dvar = [{'varname': 'delta_lim', 'varvalue':0.040}
        ]
#
# FORM method
#
column = Reliability(xvar, dvar, gfunction, None, None)
column.form(iHLRF=True, toler=1.e-6)
#
