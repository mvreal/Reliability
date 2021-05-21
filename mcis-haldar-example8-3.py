"""
Created on Mon Feb 08 17:12:00 2021
Reliabilty Analysis
Example 7.9 - Nonlinear limit state function with non-normal correlated variables
@author: MVREAL
"""
import numpy as np
from realpy import *
import time

#
# Step 0 - Column: g(Y, Z) = Y*Z-1140 = 0
#


def gfunction(x):

    g = x[0]*x[1]-1140
    return g


ti = time.time()
#
# Data input
#
# Random variables: name, probability distribution, mean and coefficient of variation

xvar = [
    {'varname': 'Y', 'vardist': 'lognormal', 'varmean': 38.00, 'varcov': 0.10, 'varhmean': 25.033776},
    {'varname': 'Z', 'vardist': 'normal', 'varmean': 54.00, 'varcov': 0.05, 'varhmean': 45.538105}
]


# Correlation matrix

corrmatrix = [[1.00, 0.30],
              [0.30, 1.00]]
#
# MCS method
#
beam = Reliability(xvar, gfunction, None, corrmatrix)
beam.mc(50_000, 1.00)
tf = time.time()
ttotal = tf - ti
print(f'Processing time = {ttotal}')
#