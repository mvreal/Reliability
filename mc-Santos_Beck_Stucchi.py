# -*- coding: utf-8 -*-
"""
Reliability of reinforced concrete beams
Santos, D. M., Stucchi, F. R., & Beck, A. T. (2014).
Reliability of beams designed in accordance with Brazilian codes.
Revista IBRACON de Estruturas e Materiais, 7(5), 723–746.
https://doi.org/10.1590/s1983-41952014000500002
"""
import numpy as np
from scipy.stats import norm
#
def gfunction(x):
    b = x[0]
    h = x[1]
    dl = x[2]
    fc = x[3]
    fy = x[4]
    g = x[5]
    q = x[7]
    thetaR = x[8]
    thetaS = x[9]
    # Função de estado limite g(x)=MR - MS = 0
    MR = thetaR * 1000. * As * fy * (
                h - dl - 0.5 * (As * fy / (0.85 * fc * b)));  # Momento de flexão resistente interno
    # O fator 1000. converte de MNm para kNm
    MS = thetaS * (g + q);  # Momento de carregamento externo (kNm)
    gx = MR - MS;  # Função estado limite
    return gx


# Parâmetro geométrico da viga
As=0.00015;    # Área de aço da seção transversal da viga (m2)

# Largura da viga: Distribuição normal
bm=0.20;         # largura da viga (m)
sb=0.0120
b=np.random.normal(bm,sb,N)
# Altura da viga: Distribuição normal
hm=0.5;        # altura da viga (m)
sh=0.0225
h=np.random.normal(hm,sh,N)
# Centroide da armadura: Distribuição lognormal
dlm=0.039
sdl=0.0110
cvdl=sdl/dlm
zetadl=np.log(1+cvdl**2)
lambdadl=np.log(dlm)-0.50*zetadl**2;
dl=np.random.lognormal(lambdadl,zetadl,N);
# Resistência a compressão do concreto: distribuição normal
fck=25; # Unidades: MPa
Vfc=0.15
fcm=1.17*fck
sfc=Vfc*fcm
fc=np.random.normal(fcm,sfc,N)

# Resistência ao escoamento do aço de reforço: distribuição lognormal
fyk=500 # Unidades: MPa
Vfy=0.05
fym=1.08*fyk
sfy=Vfy*fym
zetafy=np.log(1+Vfy**2)
lambdafy=np.log(fym)-0.50*zetafy**2;
fy=np.random.lognormal(lambdafy,Vfy,N);
#
# Momento de cálculo
#
Md=3.77 # kN.m
#
# Coeficientes de segurança
#
gamag=1.40
gamaq=1.40
k=0.10
# carga permanente:Distribuição normal
gk=Md/(gamag+gamaq*(k/(1-k))) # Unidades: kN.m
gm=gk # Unidades: kN.m
Vg=0.10
sg=Vg*gm
g=np.random.normal(gm,sg,N);

# carga acidental q: Distribuição de Gumbel 
qk=Md/(gamag*(1-k)/k+gamaq) # Unidades: kN.m
qm=0.93*qk
Vq=0.20
sq=Vq*qm
alfaq=np.sqrt((np.pi**2)/(6.*(sq)**2));
uq=qm-0.577216/alfaq;
q=np.random.gumbel(uq,1./alfaq,N);

# Erro de resistência do modelo: Distribuição lognormal
thetaRm=1.00
sthetaR=0.05
sthetaR=thetaRm*sthetaR
zetaR=np.log(1+thetaRm**2)
lambdaR=np.log(sthetaR)-0.50*zetaR**2;
thetaR=np.random.lognormal(lambdaR,zetaR,N);

# Erro de carregamento do modelo: Distribuição lognormal
thetaSm=1.00
sthetaS=0.05
sthetaS=thetaSm*sthetaS
zetaS=np.log(1+thetaSm**2)
lambdaS=np.log(sthetaS)-0.50*zetaS**2;
thetaS=np.random.lognormal(lambdaS,zetaS,N);


# Avaliação do número de falhas
index=np.where(G>=0,0,1);
Nf=sum(index);

# Avaliação do índice de confiabilidade Beta como inverso do normal 
#função de densidade cumulativa de distribuição de Pf (probabilidade de falha)
Pf=Nf/N
Beta=-norm.ppf(Pf)
print("\nBeta = {0:0.4f} e Pf = {1:0.4e}".format(Beta,Pf))