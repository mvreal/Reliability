"""
Cálculo da distribuição de probabilidade do tempo até a despassivação da armadura
Mauro Real
18/05/2023
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import norm
from scipy.special import erf
from realpy import *




#
# Step 0 - Função estado limite de despassivação: g(t) = Ccr - C(x,t)
#


def gfunction(x, d):

     #
     # Função estado limite de despassivação
     #
     t = d[0]
     D0 = d[1]
     EA = d[2]
     R = d[3]
     tl = d[4]
     t0 = x[0]
     Ccr = x[1]
     Cs = x[2]
     xc = x[3]
     Temp = x[4]
     alpha = x[5]

     # Cálculo do fator ke
     ke=np.exp(EA/R*(1./293.-1./(273.+Temp)))    
    # Cálculo do coeficiente de difussão no tempo t
     D=D0/(1.-alpha)*((1.+tl/t)**(1.-alpha)-(tl/t)**(1.-alpha))*(t0/t)**alpha*ke
    #cxtp é a concentração de cloretos em x=xc após t  anos
     Cxt = Cs * (1-erf(xc / (2 * (D*t) ** 0.5)))
     g = Ccr - Cxt

     return g


#
# Data input
#

#
# Probabilidade de falha por despassivação da armadura
#
pf = np.zeros(51)
beta = np.zeros(51)
# Tempo até a despassivação da armadura
td=np.arange(0,51)
# Dados de entrada determinísticos
D0 = 4.00*55188000.e-12 #coeficiente de difusão de cloretos em m2/anos
EA=5000.00 #EA é a ativação de energia para a difusão do cloreto [kcal/mol]
R = 1.00 #R é a constante universal dos gases perfeitos 
tl =float(28./365.) #t′ a idade do concreto quando exposto aos íons [anos]
# Geração das variáveis para as simulações de Monte Carlo
#
# Geração das variáveis aleatórias do problema
# Tempo de início da exposição t0 (anos) - distribuição normal
mediat0=28./365.
desviot0=0.036*mediat0

# Concentração crítica de cloretos - distribuição normal
mediaCcr=0.40
desvioCcr=0.10*mediaCcr

# Concentração superficial de cloretos - distribuição lognormal
mediaCs=3.50
desvioCs=0.22*mediaCs

# Cobrimento da armadura - distribuição lognormal
mediaxc=0.05
desvioxc=0.10*mediaxc

# Temperatura média anual - distribuição normal
mediaTemp=18.
desvioTemp=0.10*mediaTemp

# alpha = fator de envelhecimento do concreto - distribuição normal
mediaalpha=0.40
desvioalpha=0.10*mediaalpha

#
# Laço sobre o tempo de despassivação 
#

td[0] = 0.00
pf[0] = 0.00
beta[0] = 0.50



for i in range(1,51):
    t = td[i]
    # Random variables: name, probability distribution, mean and coefficient of variation

    xvar = [
        {'varname': 't0', 'vardist': 'lognormal', 'varmean': mediat0, 'varstd': desviot0 },
        {'varname': 'Ccr', 'vardist': 'normal', 'varmean': mediaCcr, 'varstd': desvioCcr },
        {'varname': 'Cs', 'vardist': 'lognormal', 'varmean': mediaCs, 'varstd': desvioCs },
        {'varname': 'xc', 'vardist': 'normal', 'varmean': mediaxc, 'varstd': desvioxc }, 
        {'varname': 'Temp', 'vardist': 'normal', 'varmean': mediaTemp, 'varstd': desvioTemp }, 
        {'varname': 'alpha', 'vardist': 'normal', 'varmean': mediaalpha, 'varstd': desvioalpha }
        
    ]

    # Design variables

    dvar = [
        {'varname': 't', 'varvalue': t},
        {'varname': 'D0', 'varvalue': D0},
        {'varname': 'EA', 'varvalue': EA},
        {'varname': 'R', 'varvalue': R},
        {'varname': 'tl', 'varvalue': tl}
        ]

    #
    # FORM method
    #
    vida_util = Reliability(xvar, dvar, gfunction, None)
    beta[i], xk, cos_dir, normgradyk,sigmaxneq = vida_util.form(iHLRF=True, toler=1.e-3)
    if beta[i] > 0.00:
        pf[i] = norm.cdf(-beta[i])
    else:
        pf[i] = 0.500 + norm.cdf(beta[i])
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
plt.grid()
plt.savefig('D:\Reliability\cdf_td.pdf')
plt.show()







