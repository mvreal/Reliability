# -*- coding: utf-8 -*-
"""
Reliability of reinforced concrete beams
Santos, D. M., Stucchi, F. R., & Beck, A. T. (2014).
Reliability of beams designed in accordance with Brazilian codes.
Revista IBRACON de Estruturas e Materiais, 7(5), 723–746.
https://doi.org/10.1590/s1983-41952014000500002
"""
import numpy as np
from realpy import *

#
def gfunction(x, d):
    b = x[0]
    h = x[1]
    dl = x[2]
    fc = x[3]
    fy = x[4]
    g = x[5]
    q = x[6]
    thetaR = x[7]
    thetaS = x[8]
    # Parâmetro geométrico da viga
    As1 = d[0]  # Área de aço da seção transversal da viga (m2)
    alpha_cc = d[1]
    # Função de estado limite g(x)=MR - MS = 0
    MR = thetaR * 1000. * As1 * fy * (
                h - dl - 0.5 * (As1 * fy / (alpha_cc * fc * b)))  # Momento de flexão resistente interno
    # O fator 1000. converte de MNm para kNm
    MS = thetaS * (g + q)  # Momento de carregamento externo (kNm)
    gx = MR - MS  # Função estado limite
    return gx


# Parâmetro geométrico da viga
As1 = 0.00050    # Área de aço da seção transversal da viga (m2)
#
# Momento de cálculo
#
Md = 92.00  # kN.m
#
# Coeficientes de segurança
#
gamag = 1.40
gamaq = 1.40
k = 0.2163
# carga permanente:Distribuição normal
gk = Md / (gamag + gamaq * (k / (1 - k)))  # Unidades: kN.m
gm = gk  # Unidades: kN.m
Vg = 0.10


# carga acidental q: Distribuição de Gumbel 
qk = Md / (gamag * (1 - k) / k + gamaq)  # Unidades: kN.m
qm = 0.93 * qk
Vq = 0.20


xvar = [
    {'varname': 'b', 'vardist': 'normal', 'varmean': 0.20, 'varcov': 0.06},
    {'varname': 'h', 'vardist': 'normal', 'varmean': 0.50, 'varcov': 0.045},
    {'varname': 'dl', 'vardist': 'lognormal', 'varmean': 0.039, 'varcov': 0.2821},
    {'varname': 'fc', 'vardist': 'normal', 'varmean': 1.17*25, 'varcov': 0.15},
    {'varname': 'fy', 'vardist': 'normal', 'varmean': 1.08*500, 'varcov': 0.05},
    {'varname': 'g', 'vardist': 'normal', 'varmean': gm, 'varcov': Vg},
    {'varname': 'q', 'vardist': 'gumbel', 'varmean': qm, 'varcov': Vq},
    {'varname': 'thetaR', 'vardist': 'lognormal', 'varmean': 1.00, 'varcov': 0.05},
    {'varname': 'thetaS', 'vardist': 'lognormal', 'varmean': 1.00, 'varcov': 0.05}
]

dvar = [{'varname': 'As1', 'varvalue': As1},
        {'varname': 'alpha_cc', 'varvalue': 0.85}]
#
# MCS method
#
beam = Reliability(xvar, dvar, gfunction, None, None)
beam.bucher(100, 10_000, 0.05, 1.50)
#
