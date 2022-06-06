"""
Created on Mon Feb 08 17:12:00 2021
Reliabilty Analysis
Example 7.3 - Linear limit state function with normal independent variables
@author: MVREAL
"""
from realpy import *

#
# Step 0 - Column: g(R, G, Q, W) = R-G-Q-W = 0
#


def gfunction(x, d):

    g = d[0] * x[3] * x[0] - x[1] - x[2]
    return g



# Data input
#
# Random variables: name, probability distribution, mean and coefficient of variation


xvar = [
    {'varname': 'R', 'vardist': 'normal', 'varmean': 279.81, 'varstd': 19.37},
    {'varname': 'G', 'vardist': 'normal', 'varmean': 121.50, 'varstd': 12.15},
    {'varname': 'Q', 'vardist': 'gumbel', 'varmean': 38.57, 'varstd': 9.64},
    {'varname': 'em', 'vardist': 'normal', 'varmean': 1.00, 'varcov': 0.08}
]
# Design variables

dvar = [
    {'varname': 'factor1', 'varvalue': 1.00}
   ]

#
# MC-IS method
#
column = Reliability(xvar, dvar, gfunction, None, None)
column.bucher(100, 10000, 0.05, 1.50)
#