"""
Cálculo da distribuição de probabilidade do tempo até a despassivação da armadura
Mauro Real
18/05/2023
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import norm
from scipy.special import erfinv, erfc
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
     t0 = d[4]
     Ccr = x[0]
     Cs = x[1]
     cobr = x[2]
     Temp = x[3]
     alpha = x[4]
     D0 = x[5]
     

     # Cálculo do fator ke
     ke = np.exp(EA/R*(1./293.-1./(273.+Temp)))    
     # Cálculo do coeficiente de difussão no tempo t
     D = D0/(1.-alpha)*((1.+tl/t)**(1.-alpha)-(tl/t)**(1.-alpha))*(t0/t)**alpha*ke
     # D = D0*(t0/t)**alpha*ke
     #cxtp é a concentração de cloretos em x=xc após t  anos
     xc = 2.00 * (D * t) ** 0.5  * erfinv(1. - Ccr / Cs)
     g = cobr - xc

     return g


#
# Data input
#

#
# Probabilidade de falha por despassivação da armadura
#
tf = 50
# Tempo até a despassivação da armadura
td=np.arange(0,tf+1, 5)
nt = int(tf/5)
pf = np.zeros(nt+1)
beta = np.zeros(nt+1)
data = np.zeros((nt+1,6))

# Dados de entrada determinísticos

EA=44600.00 #EA é a ativação de energia para a difusão de cloretos [J/mol]
R = 8.314 #R é a constante universal dos gases perfeitos [J/molK]
tl =float(28./365.) #t′ a idade do concreto quando exposto aos íons [anos]
t0 =float(28./365) # t0 é a idade de medida do coeficiente de difusão de cloretos

# Geração das variáveis para as simulações de Monte Carlo
#
# Geração das variáveis aleatórias do problema


# Concentração crítica de cloretos - distribuição normal
mediaCcr=0.40
desvioCcr=0.10

# Concentração superficial de cloretos - distribuição normal
mediaCs=5.50
desvioCs=1.35

# Cobrimento da armadura - distribuição normal
mediacobr=70./1000.
desviocobr=6./1000.

# Temperatura média anual - distribuição normal
mediaTemp=10.
desvioTemp=0.001

# alpha = fator de envelhecimento do concreto - distribuição normal
mediaalpha=0.40
desvioalpha=0.08

# D0 = coeficiente de difusão médio aos 28 dias = distribuição normal

mediaD0 = 6.00*31536000.e-12 #coeficiente de difusão de cloretos em m2/anos
desvioD0 = 0.64*31536000.e-12

#
# Laço sobre o tempo de despassivação 
#

td[0] = 1.00
pf[0] = 0.00
beta[0] = 100.00

i = -1

for t in td:
    i += 1
    # Random variables: name, probability distribution, mean and coefficient of variation

    xvar = [
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
        {'varname': 'tl', 'varvalue': tl},
        {'varname': 't0', 'varvalue': t0}
        ]

    #
    # FORM method
    #
    vida_util = Reliability(xvar, dvar, gfunction, None)
    beta[i], xk, cos_dir, normgradyk,sigmaxneq = vida_util.form(iHLRF=True, toler=1.e-3)
    pf[i] = norm.cdf(-beta[i])
    # Correção para quando pf>0.50, beta deve ser negativo!
    if i>=1 and pf[i-1]>0.42:
        beta[i] = -beta[i]
        pf[i] = norm.cdf(-beta[i])
    else:
        pf[0]= 0.00
        beta[0]= 6.00
    #    
    data[i,:] = xk[0:6] 
    #

# Primeiro cria um dicionário chamado res para arquivar os dados a serem inseridos no dataframe
res = {}
# Grava dados no dicionário
#
res['td'] = td
res['pf'] = pf*100.
res['Beta'] = beta
res['Ccr'] = data[:,0]
res['Cs'] = data[:,1]
res['cobr'] = data[:,2]
res['Temp'] = data[:,3]
res['alpha'] = data[:,4]
res['D0'] = data[:,5]/31536000*1e+12

# Cria então o novo dataframe, se usasse o antigo (sheet) ia gerar conflito de tamanho (número de linhas)
dfres = pd.DataFrame(res)

#===========================================
# Gravar os dados em uma planilha Excel

#-----------------------------------------------------------------------------

#   index=False evita que crie uma coluna no ínicio com o contador de linhas
#       O ExcelWriter é necessário quando mais de um dataframe é gravado no mesmo arquivo
with pd.ExcelWriter('C://Users//Mauro//OneDrive//Reliability//dados_td.xlsx', engine='openpyxl') as writer:
     dfres.to_excel(writer, sheet_name='Planilha1', index=False)    

# CDF do tempo de despassivação
# 
pf_duracon_10 = np.array([ 0.0000,  0.0000,  0.0750,  0.7470,  2.6000,  5.7310,  9.7110,  14.1450,  18.7000,  23.3480,  27.8820])
pf_duracon_20 = np.array([ 0.0000,  0.0580,  2.6530, 10.9010, 21.6100, 32.5160, 42.2720,  50.5160,  57.3410,  62.9470,  67.6820])
 
plt.plot(td,pf*100.,label="Realpy")
plt.plot(td,pf_duracon_10,label="Duracon")
plt.title('Probabilidade acumulada do tempo de despassivação')
plt.xlabel('tempo de despassivação td (anos)')
plt.ylabel('Probabilidade de falha')
plt.xlim(0,td.max())
plt.xticks(np.arange(0, max(td)+10, 10))
plt.yticks(np.arange(0, max(pf*100.)+5, 5))
plt.legend(loc='lower right', title='Gjorv - Type 1 - Temp. = 10°C')
plt.grid()
plt.savefig('D:\Reliability\cdf_td.pdf')
plt.show()






