"""
Cálculo da probabilidade de falha para uma combinação de distribuição normal e beta
Método de Monte Carlo Força Bruta
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
    {'varname': 'R', 'vardist': 'beta', 'parameter1': 0.45, 'parameter2': 1.25, 'parameter3': 0.22, 'parameter4': 0.36},
    {'varname': 'S', 'vardist': 'normal', 'varmean': 0.40, 'varstd': 0.10}      
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
column.form2(iHLRF=True, toler=1.e-6)
#






