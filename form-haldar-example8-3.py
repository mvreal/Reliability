"""
Created on Mon Feb 08 17:12:00 2021
Reliabilty Analysis
Example 7.9 - Nonlinear limit state function with non-normal correlated variables
@author: MVREAL
"""
from realpy import *

#
# Step 0 - Column: g(Y, Z) = Y*Z-1140 = 0
#


def gfunction(x):

    g = x[0]*x[1]-1140
    return g


#
# Data input
#
# Random variables: name, probability distribution, mean and coefficient of variation

xvar = [
    {'varname': 'Y', 'vardist': 'lognormal', 'varmean': 38.00, 'varcov': 0.10},
    {'varname': 'Z', 'vardist': 'normal', 'varmean': 54.00, 'varcov': 0.05}
]


# Correlation matrix

corrmatrix = [[1.00, 0.30],
              [0.30, 1.00]]
#
# FORM method
#
beam = Reliability(xvar, gfunction, None, corrmatrix)
beam.form(iHLRF=True, toler=1.e-6)
#