"""
Created on Mon Feb 08 17:12:00 2021
Reliabilty Analysis
Example 7.5 - Linear limit state function with normal correlated variables
@author: MVREAL
"""
import numpy as np
from class_reliability import *
import time

#
# Step 0 - Column: g(R, G, Q, W) = R-G-Q-W = 0
#


def gfunction(x):

    g = x[0]-x[1]-x[2]-x[3]
    return g


ti = time.time()
#
# Data input
#
# Random variables: name, probability distribution, mean and coefficient of variation


xvar = [
    {'varname': 'R', 'vardist': 'normal', 'varmean': 975.00, 'varcov': 0.15},
    {'varname': 'G', 'vardist': 'normal', 'varmean': 200.00, 'varcov': 0.07},
    {'varname': 'Q', 'vardist': 'normal', 'varmean': 300.00, 'varcov': 0.12},
    {'varname': 'w', 'vardist': 'normal', 'varmean': 150.00, 'varcov': 0.20}
]

# Correlation matrix

corrmatrix = [[1.00, 0.80, 0.00, 0.00],
              [0.80, 1.00, 0.30, 0.00],
              [0.00, 0.30, 1.00, 0.00],
              [0.00, 0.00, 0.00, 1.00]]
#
# MCS method
#
column = Reliability(xvar, gfunction, None, corrmatrix)
column.mc(1e6, 1.00)
tf = time.time()
ttotal = tf - ti
print(f'Processing time = {ttotal}')
#