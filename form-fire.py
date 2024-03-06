# -*- coding: utf-8 -*-
"""
Reliability of reinforced concrete beams in fire situation

"""
import numpy as np
from realpy import *

#
def gfunction(x, d):
    b = x[0]
    ds = x[1]
    fc = x[2]
    fy = x[3]
    Mg = x[4]
    Mq = x[5]
    thetaR = x[6]
    thetaS = x[7]
    # Parâmetro geométrico da viga
    As1 = d[0]  # Área de aço da seção transversal da viga (m2)
    alpha_cc = d[1]  # efeito Rüsch
    k_c = d[2]  # redução de fc com a temperatura
    k_s = d[3]  # redução de fy com a temperatura
    
      
    #
    # Função de estado limite g(x)=MR - MS = 0
    MR = thetaR * 1000. * As1 * fy * k_s * (
                ds - 0.5 * (As1 * fy * k_s) / (alpha_cc * fc * k_c * b))  # Momento de fletor resistente interno
    # O fator 1000. converte de MNm para kNm
    MS = thetaS * (Mg + Mq)   # Momento de carregamento externo (kNm)
    gx = MR - MS  # Função estado limite
    return gx


# Parâmetros geométricos da viga
b = 0.20 # Largura da viga em m
h = 0.50 # Altura da viga em m
As1 = 0.000942    # Área de aço da seção transversal da viga (m2)
diam_t = 5.00  # Diâmetro do estribo em mm
diam_l = 20.0  # Diâmetro da armadura longitudinal em mm
c = 30.00 # Cobrimento da armadura (até o estribo)
dl = (c + diam_t + diam_l/2.)/1000.

# Parâmetros dos materiais
fck = 25.00 # fck em MPa
alpha_cc = 1.00 # Fator do efeito Rüsch, em situação de incêncio igual a 1.00
fyk = 500.00 # fyk em MPa

# Parâmetros de incêndio
TRRF = 90.00 # Tempo requerido de resistência ao fogo
k_c = 0.6806 # Fator de redução da resistência do concreto
k_s = 0.6686 # Fator de redução da resistência do aço

#
# Momento de cálculo em situação de incêndio
#
Mrd_fi = 1000. * As1 * fyk * k_s * ( h - dl - 0.5 * (As1 * fyk * k_s) / (alpha_cc * fck * k_c * b))  # kN.m
#
# Coeficientes de segurança
#
gamag = 1.40
gamaq = 1.40
gamag_fi = 1.00
gamaq_fi = 1.00
psi_2 = 0.82
fator = 1.00
r = 2.00  # r = Mqk / Mgk
k = (gamag_fi + fator*gamaq_fi*psi_2*r)/(gamag + gamaq*r)



# carga permanente:Distribuição normal
Mgk = Mrd_fi / (k*(gamag  + gamaq*r)) # Unidades: kN.m
Mgm = 1.06 * Mgk
Vg = 0.12


# carga acidental q: Distribuição Gama
Mqk = r * Mgk  # Unidades: kN.m
Mqm = 0.21 * Mqk
Vq = 0.76


xvar = [
    {'varname': 'b', 'vardist': 'normal', 'varmean': 0.2020, 'varstd': 0.0081},
    {'varname': 'ds', 'vardist': 'normal', 'varmean': 0.4505, 'varstd': 0.0180},
    {'varname': 'fc', 'vardist': 'normal', 'varmean': 1.25*25, 'varcov': 0.17},
    {'varname': 'fy', 'vardist': 'normal', 'varmean': 1.22*500, 'varcov': 0.04},
    {'varname': 'g', 'vardist': 'normal', 'varmean': Mgm, 'varcov': Vg},
    {'varname': 'q', 'vardist': 'gamma', 'varmean': Mqm, 'varcov': Vq},
    {'varname': 'thetaR', 'vardist': 'normal', 'varmean': 1.00, 'varcov': 0.05},
    {'varname': 'thetaS', 'vardist': 'lognormal', 'varmean': 1.02, 'varcov': 0.0612}
]

dvar = [{'varname': 'As1', 'varvalue': As1},
        {'varname': 'alpha_cc', 'varvalue': 1.00},
        {'varname': 'k_c', 'varvalue': 0.6806},
        {'varname': 'k_s', 'varvalue': 0.6686}
        ]
#
# FORM method
#
beam = Reliability(xvar, dvar, gfunction, None, None)
beam.form(iHLRF=True, toler=1.e-3)
#
