# -*- coding: utf-8 -*-
"""
Reliability of reinforced concrete colums reinforced with CFRP

"""
import numpy as np
from realpy import *

#
def gfunction(x, d):
    k = d[0]
    fc = x[0]
    fy = x[1]
    Ef = x[2]
    g = x[3]
    q = x[4]

    B2 = fc
    C2 = fy
    A2 = Ef
    G6 = g
    G7 = q

    R = 230.2 - 0.476 * A2 + 21.76 * B2 + 0.1535 * C2 + (0.001073 * A2 * A2) + (0.026 * B2 * B2) + \
            (0.000077 * C2 * C2) + 0.00096 * A2 * B2 + 0.000172 * A2 * C2 + 0.00114 * B2 * C2
    S = (G6 + G7)

    gx = R - S  # Função estado limite
    return gx



xvar = [
    {'varname': 'fc', 'vardist': 'normal', 'varmean': 23.00, 'varcov': 0.15},
    {'varname': 'fyl', 'vardist': 'normal', 'varmean': 615.00, 'varcov': 0.05},
    {'varname': 'Ef', 'vardist': 'lognormal', 'varmean': 231.00, 'varcov': 0.10},
    {'varname': 'g', 'vardist': 'normal', 'varmean': 158.5594, 'varcov': 0.10},
    {'varname': 'q', 'vardist': 'gumbel', 'varmean': 317.1188, 'varcov': 0.25}
]

dvar = [{'varname': 'k', 'varvalue': 1.00}
]
#
# MCS method
#
g2r = Reliability(xvar, dvar, gfunction, None, None)
g2r.mc(100, 100_000, 0.03, 1.00)
#
