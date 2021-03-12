"""
Created on Mon Feb 08 17:12:00 2021
Reliabilty Analysis
Example 7.3 - Linear limit state function with normal independent variables
@author: MVREAL
"""
import numpy as np
from class_reliability import *
import time

#
# Step 0 - Column: g(R, G, Q, W) = R-G-Q-W = 0
#


def gfunction(x):
    a = 1.00
    b = 1.00
    c = 1.00
    d = 1.00
    g = a * x[0]- b *x[1]- c * x[2]- d * x[3]
    return g


ti = time.time()
#
# Data input
#
# Random variables: name, probability distribution, mean and coefficient of variation


xvar = [
    {'varname': 'R', 'vardist': 'normal', 'varmean': 975.00, 'varcov': 0.15},
    {'varname': 'G', 'vardist': 'lognormal', 'varmean': 200.00, 'varcov': 0.07},
    {'varname': 'Q', 'vardist': 'extvalue1', 'varmean': 300.00, 'varcov': 0.12},
    {'varname': 'w', 'vardist': 'evt3min', 'varmean': 150.00, 'varcov': 0.20}
]
#
# FORM method
#
column = Reliability(xvar, gfunction)
column.form(iHLRF=False)
tf = time.time()
ttotal = tf - ti
print(f'Processing time = {ttotal}')
#