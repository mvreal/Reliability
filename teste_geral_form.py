"""
Cálculo da probabilidade de falha para um combinação de distribuição normal e beta
Mauro Real
22/02/2023
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import norm
from scipy.stats import lognorm
from scipy.stats import beta
from scipy.special import erf
from realpy import *




#
# Step 0 - Column: g(R, S) = R - S = 0
#


def gfunction(x, d):

     #
     # Função estado limite de despassivação
     #
     a0 = d[0]
     a1 = d[1]
     r = x[0]
     s = x[1]
     g = a0 * r - a1 * s

     return g


#
# Data input
#


# Random variables: name, probability distribution, mean and coefficient of variation

xvar = [
    {'varname': 'R', 'vardist': 'normal', 'varmean': 200, 'varstd': 30 },
#    {'varname': 'S', 'vardist': 'beta', 'parameter1': 50, 'parameter2': 150, 'parameter3': 2.625, 'parameter4': 2.625}  
#     {'varname': 'R', 'vardist': 'uniform', 'parameter1': 140, 'parameter2': 260},
#     {'varname': 'S', 'vardist': 'uniform', 'parameter1': 50, 'parameter2': 150}   
#    {'varname': 'S', 'vardist': 'lognormal', 'varmean': 100, 'varstd': 20 }
    {'varname': 'S', 'vardist': 'weibull', 'varmean': 100, 'varstd': 20, 'parameter3': 50 }     
    ]

# Design variables

dvar = [
    {'varname': 'a0', 'varvalue': 1.00},
    {'varname': 'a1', 'varvalue': 1.00}
]

#
# FORM method
#
column = Reliability(xvar, dvar, gfunction, None)
column.form(iHLRF=True, toler=1.e-6)
#







