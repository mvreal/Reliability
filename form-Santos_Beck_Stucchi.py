# -*- coding: utf-8 -*-
"""
Reliability of reinforced concrete beams
Santos, D. M., Stucchi, F. R., & Beck, A. T. (2014).
Reliability of beams designed in accordance with Brazilian codes.
Revista IBRACON de Estruturas e Materiais, 7(5), 723–746.
https://doi.org/10.1590/s1983-41952014000500002
"""
import numpy as np
from class_reliability import *
import pandas as pd

#
def gfunction(x):
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
    As1 = 0.00050;  # Área de aço da seção transversal da viga (m2)
    # Função de estado limite g(x)=MR - MS = 0
    MR = thetaR * 1000. * As1 * fy * (
                h - dl - 0.5 * (As1 * fy / (0.85 * fc * b)))  # Momento de flexão resistente interno
    # O fator 1000. converte de MNm para kNm
    MS = thetaS * (g + q)  # Momento de carregamento externo (kNm)
    gx = MR - MS  # Função estado limite
    return gx
# Dimensionamento das variáveis


n = 9
xk = np.zeros(n)
alfa = np.zeros(n)
normgradyk = np.zeros(n)
sigmaxneq = np.zeros(n)

# Parâmetro geométrico da viga
As1 = 0.00050;    # Área de aço da seção transversal da viga (m2)
#
# Momento de cálculo
#
Md = 92.00 # kN.m
#
# Coeficientes de segurança
#
gamag = 1.40
gamaq = 1.40
i = -1

chi = np.linspace(0.01, 0.99, 20)
nchi = len(chi)
betak = np.zeros(nchi)

for k in chi:
    i += 1

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

    #
    # FORM method
    #
    beam = Reliability(xvar, gfunction, None, None)
    beta, xk, alpha, normgradyk, sigmaxneqk = beam.form(iHLRF=True, toler=1.e-3)
    betak[i] = beta
    #
# Primeiro cria um dicionário chamado res para arquivar os dados a serem inseridos no dataframe
res = {}
# Grava dados no dicionário
res['As'] = As1
res['gamag'] = gamag
res['gamaq'] = gamaq
res['chi'] = chi
res['Beta']= betak

# Cria então o novo dataframe, se usasse o antigo (sheet) ia gerar conflito de tamanho (número de linhas)
dfres = pd.DataFrame(res)

#===========================================
# Aqui temos dois caminhos:

#-----------------------------------------------------------------------------
# Se quiser gravar em outra aba do arquivo (pode usar o arquivo original dados.xlsx)
#   Pra usar o to_excel existe um bug na engine padrão que grava os arquivos xlsx então precisa mudar pro openpyxl
#   index=False evita que crie uma coluna no ínicio com o contador de linhas
#       O ExcelWriter é necessário quando mais de um dataframe é gravado no mesmo arquivo
with pd.ExcelWriter('Resultado.xlsx', engine='openpyxl') as writer:
#       dfres.to_excel(writer, sheet_name='Planilha2', index=False)
        dfres.to_excel(writer, sheet_name='Resultado', index=False)

