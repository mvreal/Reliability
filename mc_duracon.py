"""
Cálculo da distribuição de probabilidade do tempo até a despassivação da armadura
Mauro Real
18/05/2023
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import norm
from scipy.special import erfinv
from realpy import *




#
# Step 0 - Função estado limite de despassivação: g(t) = Ccr - C(x,t)
#


def gfunction(x, d):

     #
     # Função estado limite de despassivação
     #
     t = d[0]
     EA = d[1]
     R = d[2]
     tl = d[3]
     t0 = x[0]
     Ccr = x[1]
     Cs = x[2]
     cobr = x[3]
     Temp = x[4]
     alpha = x[5]
     D0 = x[6]

     # Cálculo do fator ke
     ke = np.exp(EA/R*(1./293.-1./(273.+Temp)))    
     # Cálculo do coeficiente de difussão no tempo t
     D = D0/(1.-alpha)*((1.+tl/t)**(1.-alpha)-(tl/t)**(1.-alpha))*(t0/t)**alpha*ke
     # D = ke*D0*(tl/t)**alpha
     #cxtp é a concentração de cloretos em x=xc após t  anos
     xc = 2.00  * erfinv(1. - Ccr / Cs)*(D*t)**0.5
     g = cobr - xc

     return g


#
# Data input
#

#
# Probabilidade de falha por despassivação da armadura
#
tf = 120
pf = np.zeros(tf+1)
beta = np.zeros(tf+1)
# Tempo até a despassivação da armadura
td=np.arange(0,tf+1)

# Dados de entrada determinísticos

EA=5000.00 #EA é a ativação de energia para a difusão do cloreto [kcal/mol]
R = 1.00 #R é a constante universal dos gases perfeitos 
tl =float(28./365.) #t′ a idade do concreto quando exposto aos íons [anos]
# Geração das variáveis para as simulações de Monte Carlo
#
# Geração das variáveis aleatórias do problema
# Tempo de início da exposição t0 (anos) - distribuição normal
mediat0=28./365.
desviot0=1./365.

# Concentração crítica de cloretos - distribuição normal
mediaCcr=0.40
desvioCcr=0.10

# Concentração superficial de cloretos - distribuição normal
mediaCs=5.50
desvioCs=1.35

# Cobrimento da armadura - distribuição normal
mediacobr=0.070
desviocobr=0.006

# Temperatura média anual - distribuição normal
mediaTemp=10.
desvioTemp=0.10

# alpha = fator de envelhecimento do concreto - distribuição normal
mediaalpha=0.40
desvioalpha=0.08

# D0 = coeficiente de difusão médio aos 28 dias = distribuição normal

mediaD0 = 6.00*31536000.e-12 #coeficiente de difusão de cloretos em m2/anos
desvioD0 = 0.64*31536000.e-12

#
# Laço sobre o tempo de despassivação 
#

td[0] = 0.00
pf[0] = 0.00
beta[0] = 0.50



for i in range(1,tf+1):
    t = td[i]
    # Random variables: name, probability distribution, mean and coefficient of variation

    xvar = [
        {'varname': 't0', 'vardist': 'normal', 'varmean': mediat0, 'varstd': desviot0 },
        {'varname': 'Ccr', 'vardist': 'normal', 'varmean': mediaCcr, 'varstd': desvioCcr },
        {'varname': 'Cs', 'vardist': 'normal', 'varmean': mediaCs, 'varstd': desvioCs },
        {'varname': 'cobr', 'vardist': 'normal', 'varmean': mediacobr, 'varstd': desviocobr }, 
        {'varname': 'Temp', 'vardist': 'normal', 'varmean': mediaTemp, 'varstd': desvioTemp }, 
        {'varname': 'alpha', 'vardist': 'normal', 'varmean': mediaalpha, 'varstd': desvioalpha },
        {'varname': 'D0', 'vardist': 'normal', 'varmean': mediaD0, 'varstd': desvioD0 }
        
    ]

    # Design variables

    dvar = [
        {'varname': 't', 'varvalue': t},
        {'varname': 'EA', 'varvalue': EA},
        {'varname': 'R', 'varvalue': R},
        {'varname': 'tl', 'varvalue': tl}
        ]

    #
    # MC method
    #
    vida_util = Reliability(xvar, dvar, gfunction, None)
    beta[i], pf[i], delta_pf, nsimul, ttotal = vida_util.mc(100, 10_000, 0.05, 1.50, igraph=False, iprint=False,)
    
    #

# Primeiro cria um dicionário chamado res para arquivar os dados a serem inseridos no dataframe
res = {}
# Grava dados no dicionário
#
res['td'] = td
res['pf'] = pf
res['beta'] = beta

# Cria então o novo dataframe, se usasse o antigo (sheet) ia gerar conflito de tamanho (número de linhas)
dfres = pd.DataFrame(res)

#===========================================
# Gravar os dados em uma planilha Excel

#-----------------------------------------------------------------------------

#   index=False evita que crie uma coluna no ínicio com o contador de linhas
#       O ExcelWriter é necessário quando mais de um dataframe é gravado no mesmo arquivo
with pd.ExcelWriter('D:\Reliability\dados_td.xlsx', engine='openpyxl') as writer:
     dfres.to_excel(writer, sheet_name='Planilha1', index=False)    

# Histograma do tempo de despassivação
# 
plt.plot(td,pf)
plt.title('Probabilidade acumulada do tempo de despassivação')
plt.xlabel('tempo de despassivação td (anos)')
plt.ylabel('Probabilidade de falha')
plt.xlim(0,td.max())
plt.ylim(0,1.00)
plt.grid()
plt.savefig('D:\Reliability\cdf_td.pdf')
plt.show()







