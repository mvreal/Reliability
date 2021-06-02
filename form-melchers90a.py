"""
Created on Tue Apr 12 14:29:00 2021
Reliabilty Analysis
Example Melchers - 1990
MELCHERS, R.E. Search based importance sampling.
Structural Safety, 9(1990) 117-128
@author: MVREAL

Equation 7a

"""
from realpy import *

#
# Step 0 - Limit state funcion
#

def gfunction(x, d):

    g = d[0] * x[0] - d[1] * x[1] - d[2] * x[2] - d[3] * x[3] + d[4]
    return g

#
# Data input
#
# Random variables: name, probability distribution, mean and coefficient of variation

xvar = [
    {'varname': 'X1', 'vardist': 'normal', 'varmean': 0.00, 'varstd': 1.00},
    {'varname': 'X2', 'vardist': 'normal', 'varmean': 0.00, 'varstd': 1.00},
    {'varname': 'X3', 'vardist': 'normal', 'varmean': 0.00, 'varstd': 1.00},
    {'varname': 'X4', 'vardist': 'normal', 'varmean': 0.00, 'varstd': 1.00},
]


dvar = [
    {'varname': 'factor1', 'varvalue': 1.00},
    {'varname': 'factor2', 'varvalue': 1.00},
    {'varname': 'factor3', 'varvalue': 1.00},
    {'varname': 'factor4', 'varvalue': 1.00},
    {'varname': 'constant', 'varvalue': 6.00}

]
#
# FORM method
#
test = Reliability(xvar, dvar, gfunction)
test.form(iHLRF=True, toler=1.e-6)
#