"""
Created on Mon Feb 08 17:12:00 2021
Reliabilty Analysis
Example 7.7 - Linear limit state function with non-normal independent variables
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
    {'varname': 'R', 'vardist': 'lognormal', 'varmean': 975.00, 'varcov': 0.15},
    {'varname': 'G', 'vardist': 'normal', 'varmean': 200.00, 'varcov': 0.07},
    {'varname': 'Q', 'vardist': 'gumbel', 'varmean': 300.00, 'varcov': 0.12},
    {'varname': 'w', 'vardist': 'gumbel', 'varmean': 150.00, 'varcov': 0.20}
]
#
# MC-IS adaptative method
#
column = Reliability(xvar, gfunction)
column.bucher(100, 5000, 0.03, 1.5)
#