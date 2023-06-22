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



def gfunction(x, d, t):

     #
     # Função estado limite de despassivação
     #
     
     EA = d[0]
     R = d[1]
     tl = d[2]
     t0 = d[3]
     Ccr = x[0]
     Cs = x[1]
     cobr = x[2]
     Temp = x[3]
     alpha = x[4]
     D0 = x[5]
     fator_D = x[6]
     

     # Cálculo do fator ke
     ke = np.exp(EA/R*(1./293.-1./(273.+Temp)))   
     # Cálculo do coeficiente de difussão no tempo tlim
     #tlim = 30.00 # tempo limite em anos
     #Dlim = fator_D*D0/(1.-alpha)*((1.+tl/tlim)**(1.-alpha)-(tl/tlim)**(1.-alpha))*(t0/tlim)**alpha*ke 
     # Cálculo do coeficiente de difussão no tempo t
     D = fator_D*D0/(1.-alpha)*((1.+tl/t)**(1.-alpha)-(tl/t)**(1.-alpha))*(t0/t)**alpha*ke
     #
     #if D<Dlim:
     #     D = Dlim
     #cxtp é a concentração de cloretos em x=xc após t  anos
     xc = 2.00 * (D * t) ** 0.5  * erfinv(1. - Ccr / Cs)
     g = cobr - xc

     return g, D, xc


#
# Data input
#

x = []
d = []
#
# Probabilidade de falha por despassivação da armadura
#

# Dados de entrada determinísticos

EA=35000.00 #EA é a ativação de energia para a difusão de cloretos [J/mol]
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
mediaTemp=20.
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


# Random variables: name, probability distribution, mean and coefficient of variation

xvar = [
    {'varname': 'Ccr', 'vardist': 'normal', 'varmean': mediaCcr, 'varstd': desvioCcr },
    {'varname': 'Cs', 'vardist': 'normal', 'varmean': mediaCs, 'varstd': desvioCs },
    {'varname': 'cobr', 'vardist': 'normal', 'varmean': mediacobr, 'varstd': desviocobr }, 
    {'varname': 'Temp', 'vardist': 'normal', 'varmean': mediaTemp, 'varstd': desvioTemp }, 
    {'varname': 'alpha', 'vardist': 'normal', 'varmean': mediaalpha, 'varstd': desvioalpha },
    {'varname': 'D0', 'vardist': 'normal', 'varmean': mediaD0, 'varstd': desvioD0 },
    {'varname': 'fator_D', 'vardist': 'normal', 'varmean': 0.832, 'varstd': 0.024 }
]

# Design variables
dvar = [
    {'varname': 'EA', 'varvalue': EA},
    {'varname': 'R', 'varvalue': R},
    {'varname': 'tl', 'varvalue': tl},
    {'varname': 't0', 'varvalue': t0}
    ]

#
# Generation of the random variables
#
ns = 10_000
nxvar = len(xvar)
ndvar = len(dvar)

Df = np.zeros(ns)
xcr = np.zeros(ns)
vida_util = Reliability(xvar, dvar, gfunction, None)
xmatrix = vida_util.generator(ns)
#    
#
# Monte Carlo simulations
#
#
# Probabilidade de falha por despassivação da armadura
#
tf = 50
# Tempo até a despassivação da armadura
td=np.arange(0,tf+1, 5)
nt = int(tf/5)
pf = np.zeros(nt+1)
td[0] = 1.00
pf[0] = 0.00
igx = np.zeros(ns)
k = -1
for t in td:
    k += 1

    for i in range(ns):
        
        x[0:nxvar] = xmatrix[i,0:nxvar]
        d = [item['varvalue'] for item in dvar[:4]]
        gx, Dt, xc = gfunction(x,d,t)
        igx[i] = np.where(gx <= 0.00, 1, 0)
          
    
    nfail = sum(igx)
    pf[k] = nfail/ns

# Impressão dos resultados
# 


# Primeiro cria um dicionário chamado res para arquivar os dados a serem inseridos no dataframe
res = {}
# Grava dados no dicionário
#

res['Ccr'] = xmatrix[:,0] # %
res['Cs'] = xmatrix[:,1] # %
res['cobr'] = xmatrix[:,2]*1000 # mm
res['Temp'] = xmatrix[:,3] # graus Celsius
res['alpha'] = xmatrix[:,4] # sem unidade
res['D0'] = xmatrix[:,5]/31.536e6 # m2/s
res['fator'] = xmatrix[:,6] # sem unidade
res['Dt'] = Df/31.536e6 # m2/s
res['xcr'] = xcr*1000 # mm



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
plt.plot(td,pf_duracon_20,label="Duracon")
plt.title('Probabilidade acumulada do tempo de despassivação')
plt.xlabel('tempo de despassivação td (anos)')
plt.ylabel('Probabilidade de falha')
plt.xlim(0,td.max())
plt.xticks(np.arange(0, max(td)+10, 10))
plt.yticks(np.arange(0, max(pf*100.)+5, 5))
plt.legend(loc='lower right', title='Gjorv - Type 1 - Temp. = 20°C')
plt.grid()
plt.savefig('C://Users//Mauro//OneDrive//Reliability//cdf_td.pdf')
plt.show()

"""       
#
# Histograma de t0
# 
hx, hy, _  = plt.hist(t0*365, bins=100, density=True,color="blue")
plt.title('Histograma do tempo de início')
plt.xlabel('tempo de início t0 (dias)')
plt.ylabel('frequência')
plt.grid()
plt.savefig('./histogram_t0.pdf')
plt.show()
#
# Histograma da temperatura
# 
hx, hy, _  = plt.hist(Temp, bins=100, density=True,color="blue")
plt.title('Histograma da temperatura')
plt.xlabel('Temperatura (ºC)')
plt.ylabel('frequência')
plt.grid()
plt.savefig('./histogram_Temp.pdf')
plt.show()
"""
Ccr = xmatrix[:,0] # %
Cs = xmatrix[:,1] # %
cobr = xmatrix[:,2]*1000 # mm
Temp= xmatrix[:,3] # graus Celsius
alpha= xmatrix[:,4] # sem unidade
D0 = xmatrix[:,5]/31.536e6 # m2/s
#
# Histograma da concentração superficial de cloretos
# 
hx, hy, _  = plt.hist(Cs, bins=100, density=True,color="blue")
plt.title('Histograma da concentração superficial Cs')
plt.xlabel('Concentração superficial de cloretos (%)')
plt.ylabel('frequência')
plt.grid()
plt.savefig('./histogram_cs.pdf')
plt.show()
#
# Histograma da concentração crítica de cloretos
# 
hx, hy, _  = plt.hist(Ccr, bins=100, density=True,color="blue")
plt.title('Histograma concentração crítica de cloretos (Ccr)')
plt.xlabel('concentração crítica de cloretos (%)')
plt.ylabel('frequência')
plt.grid()
plt.savefig('./histogram_ccr.pdf')
plt.show()
#
# Histograma do cobrimento da armadura
# 
hx, hy, _  = plt.hist(cobr, bins=100, density=True,color="blue")
plt.title('Histograma do cobrimento da armadura (xc)')
plt.xlabel('cobrimento da armadura (mm)')
plt.ylabel('frequência')
plt.grid()
plt.savefig('./histogram_xc.pdf')
plt.show()
#
# Histograma do fator de envelhecimento do concreto
# 
hx, hy, _  = plt.hist(alpha, bins=100, density=True,color="blue")
plt.title(r'Histograma do fator de envelhecimento $\alpha$')
plt.xlabel(r'fator de envelhecimento $\alpha$')
plt.ylabel('frequência')
plt.grid()
plt.savefig('./histogram_alpha.pdf')
plt.show()
#
# Histograma do coeficiente de difusão no tempo de despassivação
# 
hx, hy, _  = plt.hist(D0, bins=100, density=True,color="blue")
plt.title('Histograma do coeficiente de difusão inicial')
plt.xlabel('Coeficiente de difusão final (m2/s)')
plt.ylabel('frequência')
plt.grid()
plt.savefig('./histogram_D.pdf')
plt.show()
