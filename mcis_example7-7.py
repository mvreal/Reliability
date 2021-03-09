"""
Created on Mon Feb 08 17:12:00 2021
Reliabilty Analysis
Example 7.7 - Linear limit state function with non-normal independent variables
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
    {'varname': 'R', 'vardist': 'lognormal', 'varmean': 975.00, 'varcov': 0.15, 'varhmean': 699.675868},
    {'varname': 'G', 'vardist': 'normal', 'varmean': 200.00, 'varcov': 0.07, 'varhmean': 204.067350},
    {'varname': 'Q', 'vardist': 'gumbel', 'varmean': 300.00, 'varcov': 0.12, 'varhmean': 328.647738},
    {'varname': 'w', 'vardist': 'gumbel', 'varmean': 150.00, 'varcov': 0.20, 'varhmean': 166.960779}
]
#
# MCS method
#
column = Reliability(xvar, gfunction)
column.mc(50_000, 1.00)
tf = time.time()
ttotal = tf - ti
print(f'Processing time = {ttotal}')
#