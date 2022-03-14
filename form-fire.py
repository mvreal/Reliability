# -*- coding: utf-8 -*-
"""
Reliability of reinforced concrete beams in fire situation

"""
import numpy as np
from realpy import *

#
def gfunction(x, d):
    b = x[0]
    h = x[1]
    c = x[2]
    fc = x[3]
    fy = x[4]
    g = x[5]
    q = x[6]
    thetaR = x[7]
    thetaS = x[8]
    # Parâmetro geométrico da viga
    As1 = d[0]  # Área de aço da seção transversal da viga (m2)
    alpha_cc = d[1]  # efeito Rüsch
    alpha_c = d[2]  # redução de fc com a temperatura
    alpha_y = d[3]  # redução de fy com a temperatura
    diam_l = d[4]
    L = d[5]
    #
    dl = c + diam_l/2000.
    #
    # Função de estado limite g(x)=MR - MS = 0
    MR = thetaR * 1000. * As1 * fy * alpha_y * (
                h - dl - 0.5 * (As1 * fy * alpha_y) / (alpha_cc * fc * alpha_c * b))  # Momento de fletor resistente interno
    # O fator 1000. converte de MNm para kNm
    MS = thetaS * (g + q) * L ** 2 / 8.0  # Momento de carregamento externo (kNm)
    gx = MR - MS  # Função estado limite
    return gx


# Parâmetro geométrico da viga
As1 = 0.00050    # Área de aço da seção transversal da viga (m2)
diam_t = 0.00  # Diâmetro do estribo em mm
diam_l = 12.5  # Diâmetro da armadura longitudinal em mm
#
# Momento de cálculo
#
Msd = 92.00  # kN.m
#
# Coeficientes de segurança
#
gamag = 1.40
gamaq = 1.40
r = 2.00  # r = gm / qm
L = 5.00  # vão L = 5 m

# carga acidental q: Distribuição de Gumbel
qm = 8.00 / L ** 2 * Msd / (r * gamag / 1.05 + gamaq)  # Unidades: kN/m
qmi = 0.24 * qm
Vq = 0.65


# carga permanente:Distribuição normal
gm = r * qm  # Unidades: kN.m
gmi = 1.05 * gm
Vg = 0.10


xvar = [
    {'varname': 'b', 'vardist': 'normal', 'varmean': 0.20, 'varstd': 0.012},
    {'varname': 'h', 'vardist': 'normal', 'varmean': 0.50, 'varstd': 0.0225},
    {'varname': 'c', 'vardist': 'lognormal', 'varmean': 0.035, 'varstd': 0.011},
    {'varname': 'fc', 'vardist': 'normal', 'varmean': 1.17*25, 'varcov': 0.10},
    {'varname': 'fy', 'vardist': 'normal', 'varmean': 1.08*500, 'varcov': 0.08},
    {'varname': 'g', 'vardist': 'normal', 'varmean': gmi, 'varcov': Vg},
    {'varname': 'q', 'vardist': 'gumbel', 'varmean': qmi, 'varcov': Vq},
    {'varname': 'thetaR', 'vardist': 'lognormal', 'varmean': 1.00, 'varcov': 0.05},
    {'varname': 'thetaS', 'vardist': 'lognormal', 'varmean': 1.00, 'varcov': 0.05}
]

dvar = [{'varname': 'As1', 'varvalue': As1},
        {'varname': 'alpha_cc', 'varvalue': 0.85},
        {'varname': 'alpha_c', 'varvalue': 0.60},
        {'varname': 'alpha_y', 'varvalue': 0.78},
        {'varname': 'diam_l', 'varvalue': diam_l},
        {'varname': 'L', 'varvalue': L},
        ]
#
# FORM method
#
beam = Reliability(xvar, dvar, gfunction, None, None)
beam.form(iHLRF=True, toler=1.e-6)
#
