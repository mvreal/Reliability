"""
Created on Mon Feb 08 17:12:00 2021
Reliabilty Analysis
@author: MVREAL
"""
from realpy import *

#
# Step 0 - Bridge bending ultimate limit state: g(R, D1, D2, D3,Mtt) = R-D1-D2-D3-Mtt = 0
#


def gfunction(x, d):

    g = d[0] * x[0] - d[1] * x[1] - d[2] * x[2] - d[3] * x[3] - d[4] * x[4]
    return g


#
# Data input
#
# Random variables: name, probability distribution, mean and coefficient of variation


xvar = [
    {'varname': 'R', 'vardist': 'lognormal', 'varmean': 27038, 'varcov': 2529/27038},
    {'varname': 'D1', 'vardist': 'normal', 'varmean': 6171, 'varcov': 494/6171},
    {'varname': 'D2', 'vardist': 'normal', 'varmean': 542, 'varcov': 54/542},
    {'varname': 'D3', 'vardist': 'normal', 'varmean': 545, 'varcov': 136/545},
    {'varname': 'Mtt', 'vardist': 'normal', 'varmean': 4527, 'varcov': 888/4527}
]
# Design variables

dvar = [
    {'varname': 'factor1', 'varvalue': 1.00},
    {'varname': 'factor2', 'varvalue': 1.00},
    {'varname': 'factor3', 'varvalue': 1.00},
    {'varname': 'factor4', 'varvalue': 1.00},
    {'varname': 'factor5', 'varvalue': 1.00}
]
#
# FORM method
#
bridge = Reliability(xvar, dvar, gfunction)
bridge.form(iHLRF=True)

#