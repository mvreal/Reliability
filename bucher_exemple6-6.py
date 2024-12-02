"""
Cálculo da probabilidade de falha de um grupo de estacas

Created on Mon Dec 02 10:30:00 2024
Reliabilty Analysis
Example 6-5 - Reliability Analysis of a Truss Bar
@author: MVREAL
"""
from realpy import *


def gfunction(x, d):

    fy = x[0]
    G = x[1]
    Q = x[2]
    Ast = d[0]
    fator_kN = d[1]

    g = fator_kN*Ast*fy - (G + Q)
    return g

#
# Data input
#

# Random variables: name, probability distribution, mean and coefficient of variation


xvar = [
    {'varname': 'fy', 'vardist': 'normal', 'varmean': 275.00, 'varstd': 0.10*275.00},
    {'varname': 'G', 'vardist': 'normal', 'varmean': 60.00, 'varstd': 0.10*60},
    {'varname': 'Q', 'vardist': 'normal', 'varmean': 70.00, 'varstd': 0.30*70},
]

# Design variables

dvar = [
    {'varname': 'Ast', 'varvalue': 804.25},
    {'varname': 'fator_kN', 'varvalue': 1./1000.00},
]

#
# MCS method
#
platform = Reliability(xvar, dvar, gfunction, None, None)
platform.bucher2(100, 10000, 0.05, 1.50)
#