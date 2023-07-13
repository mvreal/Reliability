"""
Cálculo da probabilidade de falha para uma combinação de distribuição normal e beta
Método de Monte Carlo Força Bruta
22/02/2023
"""
import numpy as np
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
    {'varname': 'R', 'vardist': 'normal', 'varmean': 25, 'varstd': 0.01 },
    {'varname': 'S', 'vardist': 'gamma', 'varmean': 15, 'varcov': 0.25 }
    ]

# Design variables

dvar = [
    {'varname': 'a0', 'varvalue': 1.00},
    {'varname': 'a1', 'varvalue': 1.00}
]

#
# MCS  method
#
column = Reliability(xvar, dvar, gfunction)
column.mc(10, 100000, 0.05, 1.00)
#







