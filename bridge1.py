"""
Created on Mon Feb 08 17:12:00 2021
Reliabilty Analysis
@author: MVREAL
"""
import numpy as np
from realpy import *
import time

#
# Step 0 - Bridge bending ultimate limit state: g(R, D1, D2, D3,Mtt) = R-D1-D2-D3-Mtt = 0
#


def gfunction(x):

    g = x[0] - x[1] - x[2] - x[3] - x[4]
    return g


ti = time.time()
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
#
# FORM method
#
bridge = Reliability(xvar, gfunction)
bridge.form(iHLRF=True)
tf = time.time()
ttotal = tf - ti
print(f'Processing time = {ttotal}')
#