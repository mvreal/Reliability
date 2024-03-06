"""
Created on Thu Dec 08 15:44:00 2023
Reliabilty Analysis
Example of probabilistic column stability analysis 
Steel column with an I cross-section
@author: MVREAL
"""
from realpy import *
import math

#
# Step 0 - Beam: g(Y, Z, M) = Y*Z-M = 0
#


def gfunction(x, d):
    k = d[0]
    b = d[1]
    h = d[2]
    le = d[3]
    E = x[0]
    tw = x[1]
    tf = x[2]
    P = x[3]
    pi = math.pi
    I = tw*h**3/12. + 2.*(b*tf**3/12. + (h/2.-tf/2.)**2*(b*tf))
    Pcr = pi**2*E*I/(k*le)**2
    g = Pcr - P
    return g


#
# Data input
#
# Random variables: name, probability distribution, mean and coefficient of variation


xvar = [
    {'varname': 'E', 'vardist': 'lognormal', 'varmean': 200.00e3, 'varcov': 0.10},
    {'varname': 'tw', 'vardist': 'normal', 'varmean': 5.00, 'varcov': 0.05},
    {'varname': 'tf', 'vardist': 'normal', 'varmean': 5.00, 'varcov': 0.05},
    {'varname': 'P', 'vardist': 'gumbel', 'varmean': 500.00e3, 'varcov': 0.20}
]

# Design variables

dvar = [
    {'varname': 'k', 'varvalue': 1.00},
    {'varname': 'b', 'varvalue': 100.00},
    {'varname': 'h', 'varvalue': 200.00},
    {'varname': 'le', 'varvalue': 5000.00}
]

# MC-IS adaptative method
#
beam = Reliability(xvar, dvar, gfunction)
beam.mc2(1, 1000000, 0.05, 1.00)
#