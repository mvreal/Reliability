"""
Created on Mon Feb 08 17:12:00 2021
Reliabilty Analysis
Example 7.3 - Linear limit state function with normal independent variables
@author: MVREAL
"""
from class_reliability import *

#
# Step 0 - Column: g(R, G, Q, W) = R-G-Q-W = 0
#


def gfunction(x):

    g = x[0]-x[1]-x[2]-x[3]
    return g


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
#
# MCS method
#
column = Reliability(xvar, gfunction)
column.mc(50, 10000, 0.01, 1.00)
#