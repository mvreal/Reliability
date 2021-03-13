"""
Created on Mon Feb 08 17:12:00 2021
Reliabilty Analysis
Example 7.9 - Nonlinear limit state function with non-normal correlated variables
@author: MVREAL
"""
import numpy as np
from class_reliability import *
import time

#
# Step 0 - Column: g(Y, Z, M) = Y*Z-M = 0
#


def gfunction(x):

    g = x[0]*x[1]-x[2]
    return g


ti = time.time()
#
# Data input
#
# Random variables: name, probability distribution, mean and coefficient of variation

xvar = [
    {'varname': 'Y', 'vardist': 'lognormal', 'varmean': 40.00, 'varcov': 0.125},
    {'varname': 'Z', 'vardist': 'weibull', 'varmean': 50.00, 'varcov': 0.05,'varinf':0.00},
    {'varname': 'M', 'vardist': 'gumbel', 'varmean': 1000.00, 'varcov': 0.20}
]


# Correlation matrix

corrmatrix = [[1.00, 0.40, 0.00],
              [0.40, 1.00, 0.00],
              [0.00, 0.00, 1.00]]
#
# FORM method
#
beam = Reliability(xvar, gfunction, None, corrmatrix)
beam.form(iHLRF=True)
tf = time.time()
ttotal = tf - ti
print(f'Processing time = {ttotal}')
#