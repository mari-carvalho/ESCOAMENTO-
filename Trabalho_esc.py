#------------------------Fase óleo------------------------#
import np as np
import numpy as np
import math as mt

# Ponto de Bolha
do = rho_o/rho_w
API = 141.5/do - 131.5

dg = rho_gas/rho_ar   # em 101.325kPa e 15Celsius Par = 1.225kg/m³
#dgn = dg*(1 + (5.912*(10**-5))*API*Tsep*np.log(Psep/114.7))       # Tsep em F

a_pb = (7.916*10**-4)*(API**1.5410)-(4.561*10**-5)*(T**1.3911)   # T em F

# Razão de Solubilidade gás-óleo
# P <= Pb
a_rs = (7.916*10**-4)*(API**1.541)+(-4.561*10**-5)*(T**1.3911)   # T em F
rs = (((P/112.727)+12.34)*(dg**0.8439)*(10**a_rs))**1.73184   # P em psia e T em F

pb = ((112.727*(rs**0.577421))/((dg**0.8439)*10**a_pb)) - 1391.051   # T em F

# Fator Volume-Formação
# P < Pb
bo = 1.0113+(7.2046*10**-5)*((rs**0.3738)*(dg**0.2914/do**0.6265)+0.24626*(T**0.5371))**3.0936)   # T em F

# Massa Específica do Óleo

# Compressibilidade Isotérmica do Óleo
# P >= Pb
co = 1.705*(10**-7)*(rs**0.69357)*(dg**0.1885)*(API**0.3272)*(T**0.6729)*(P**-0.5906)   # T em F

# Viscosidade do Óleo Morto
a_uod = 10**(0.43+8.33/API)

uod = 0.32+(1.8*10**7/API**4.53)*((360/T-260)**a_uod)   # T em Rankine

# Viscosidade do Óleo Saturado
# P <= Pb
a_uob = 10**(-7.4*10**-4*rs+2.2*10**-7*rs**2)
b_uob = 0.68/(10**8.62*10**-5*rs) + 0.25/(10**1.1*10**-3*rs) + 0.062/(10**3.74*10**-3*rs)

uob = a_uob*uod**b_uob

# Viscosidade do Óleo Sub-Saturado
uo = uob+(0.001*(P-Pb))*(0.024*(uob**1.6)+0.038*(uob**0.56))


#------------------------Fase Gás------------------------#

# Fator de Compressibilidade Isotérmico
a_z = 1.39*((Tpr-0.92)**1/2)-0.36*Tpr-0.101
b_z = (0.62-0.23*Tpr)*Ppr+((0.066/(Tpr-0.86))-0.037)*Ppr**2+(0.32*Ppr**6/(10**(9*(Tpr-1))))
c_z = 0.132-0.32*np.log(Tpr)    # log 10????
d_z = np.antilog(0.3106-0.49*Tpr+0.1824*Tpr**2)   # antilog???

z = a_z+((1-a_z)/(2.71**b_z))+c_z*Ppr**d_z   # euler ???

# Massa Específica do Gás
rho_g = P*Mg/(z*R*T)    # T em K

# Viscosidade do Gás
x_ug = 3.47+(1588/T)+0.0009*Mg   # T em Rankine
y_ug = 1.66378-0.04679*x_ug
a_ug = x_ug*rho_g**y_ug
k_ug = (0.807*(Tpr**0.618)-0.357*np.exp(-0.449*Tpr)+0.34*np.exp(-4.058*Tpr)+0.018)/(0.9490*(Tpc/((Mg**3)*(Ppc**4)))**1/6))   # T em Rankine

ug = 10**-4*k_ug*np.exp(a_ug)

# Fator Volume Formação
bg = Psc/Tsc*(z*T/P)   # Psc = 14.7 psia e Tsc = 60 F

# Compressibilidade Isotérmica
cg = 1/P-1/Z*(dz/dp)


#------------------------Fase Água------------------------#
