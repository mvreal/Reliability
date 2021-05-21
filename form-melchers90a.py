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


def gfunction(x):

    g = x[0] + x[1] - x[2] - x[3] + 6.00
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
#
# FORM method
#
test = Reliability(xvar, gfunction)
test.form(iHLRF=True, toler=1.e-6)
#