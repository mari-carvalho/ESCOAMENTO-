# Fase óleo

import numpy as np
import math as mt

# Ponto de Bolha
do = rho_o/rho_w
API = 141.5/do - 131.5
if API <= 30:
    c1_pb = 27.624
    c2_pb = 0.914328
    c3_pb = 11.172
elif API > 30:
        c1_pb = 56.18
        c2_pb = 0.84246
        c3_pb = 10.393

dg = rho_gas/rho_ar      # em 101.325kPa e 15Celsius Par = 1.225kg/m³
dgn = dg*(1 + (5.912*(10**-5))*API*Tsep*np.log(Psep/114.7))       # Tsep em F
a = -c3_pb*API/T       # temperatura em Rankine

# Razão de Solubilidade gás-óleo
# P < Pb
if API <= 30:
    c1_rs = 0.0362
    c2_rs = 1.0937
    c2_rs = 25.7240
elif API > 30:
    c1_rs = 0.0178
    c2_rs = 1.1870
    c3_rs = 23.931

rs = (c1_rs*dgn*(P)**c2_rs)*np.exp(c2_rs*API/T)   # P em psia e T em Rankine

pb = ((c1_pb*rs/dgn)*10**a)**c2_pb

# Fator Volume-Formação
# P < Pb
if API <= 30:
    c1_bo = 4.677*10**-4
    c2_bo = 1.751*10**-5
    c2_bo = -1.811*10**-8
elif API > 30:
    c1_bo = 4.670*10**-4
    c2_bo = 1.100*10**-5
    c3_bo = 1.337*10**-9

bo = 1 + c1_bo*rs + (T-520)*(API/dgn)*(c2_bo+c3_bo*rs)   # T em Rankine

# Massa Específica do Óleo

# Compressibilidade Isotérmica do Óleo
# P >= Pb
co = (-1433 + 5*rsb + 17.2*T - 1180*dgn + 12.61*API)/(P*10**5)   # T em F

# Viscosidade do Óleo Morto
a_uod = 10**(3.0324-0.02023*API)

uod = (10**(a_uod*(T)**-1.163)) -1   # T em F

# Viscosidade do Óleo Saturado
# P <= Pb
a_uob = 10.715*(rs+100)**-0.515
b_uob = 5.44*(rs+150)**-0.338

uob = a*uod**b

# Viscosidade do Óleo Sub-Saturado
m_uo = 2.6*P**1.187*np.exp(-11.513-(8.98*10**-5)*P)

uo = uob*(P/Pb)**m_uo


# Fase Gás

# Fator de Compressibilidade Isotérmico
a_z = 1.39*((Tpr-0.92)**1/2)-0.36*Tpr-0.101
b_z = (0.62-0.23*Tpr)*Ppr+((0.066/(Tpr-0.86))-0.037)*Ppr**2+(0.32*Ppr**6/(10**(9*(Tpr-1))))
c_z = 0.132-0.32*np.log(Tpr)    # log 10????
d_z = np.antilog(0.3106-0.49*Tpr+0.1824*Tpr**2)   # antilog???

z = a_z+((1-a_z)/(2.71**b_z))+c_z*Ppr**d_z   # euler ???

# Massa Específica do Gás
rho_g = P*Mg/(z*R*T)    # T em K

