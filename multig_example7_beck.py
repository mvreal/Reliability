"""
Created on Tue Apr 15 11:13:00 2021
Reliabilty Analysis
Example Beck - 2019
BECK, A.T. Confiabilidade e Segurança das Estruturas.
Elsevier, ISBN 978-85-352-8895-7.
@author: MVREAL
"""
from realpy import *

#
# Step 0 - Limit state funcion
#


def gR1(x):
    r1 = 4.00
    a1 = np.pi * r1 ** 2
    v = 300
    k = 1 / 2
    h = k * v
    le = np.sqrt(v ** 2 + h ** 2)

    g = a1*x[0] / 1000 - le / (2 * v) * (x[2] - x[3] / k)

    return g

def gE1(x):
    r1 = 4.00
    i1 = (np.pi * r1 ** 4) / 4
    v = 300
    k = 1 / 2
    h = k * v
    le = np.sqrt(v ** 2 + h ** 2)

    g = np.pi ** 2 * x[1] * i1 / le ** 2 - le / (2 * v) * (-x[2] + x[3] / k)

    return g


def gE2(x):
    r2 = 5.20
    i2 = (np.pi * r2 ** 4) / 4
    v = 300
    k = 1 / 2
    h = k * v
    le = np.sqrt(v ** 2 + h ** 2)

    g = np.pi ** 2 * x[1] * i2 / le ** 2 - le / (2 * v) * (+x[2] + x[3] / k)

    return g

#
# Data input
#

# Random variables: name, probability distribution, mean and coefficient of variation


xvar = [
    {'varname': 'S', 'vardist': 'normal', 'varmean': 24.5643, 'varcov': 0.10},
    {'varname': 'E', 'vardist': 'normal', 'varmean': 70.00, 'varcov': 0.03},
    {'varname': 'H', 'vardist': 'normal', 'varmean': 2.00, 'varcov': 0.20},
    {'varname': 'V', 'vardist': 'normal', 'varmean': 1.00, 'varcov': 0.20},
]
glist = [gR1, gE1, gE2]
test = Reliability(xvar, glist, None, None)
test.multig(xvar, glist)