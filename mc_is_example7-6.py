"""
Created on Mon Feb 08 17:12:00 2021
Reliabilty Analysis
Example 7.6 - Nonlinear limit state function with normal correlated variables
@author: MVREAL
"""
import numpy as np
from class_reliability import *
import time

#
# Step 0 - Beam: g(Y, Z, M) = Y*Z-M = 0
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
    {'varname': 'Y', 'vardist': 'normal', 'varmean': 40.00, 'varcov': 0.125, 'varhmean': 28.877084},
    {'varname': 'Z', 'vardist': 'normal', 'varmean': 50.00, 'varcov': 0.05, 'varhmean': 46.471156},
    {'varname': 'M', 'vardist': 'normal', 'varmean': 1000.00, 'varcov': 0.20, 'varhmean': 1341.995396}
]


# Correlation matrix

corrmatrix = [[1.00, 0.40, 0.00],
              [0.40, 1.00, 0.00],
              [0.00, 0.00, 1.00]]
#
# MCS method
#
beam = Reliability(xvar, gfunction, None, corrmatrix)
beam.mc(50_000, 1.00)
tf = time.time()
ttotal = tf - ti
print(f'Processing time = {ttotal}')
#