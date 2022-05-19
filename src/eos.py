from constants import *

def eos_rhoT(rho,T,Xi):

    #composition
    X = Xi['H']
    Y = Xi['He']
    Z = Xi['Z']

    #mean molecular weight
    mu_inv = X + 0.25*Y + 0.0625*Z
    mu = 1./mu_inv

    #ideal gas
    Pgas = (rho/mu)*cgas*T
    Egas = 1.5*Pgas/rho
    Sgas = (Pgas/rho + Egas)/T

    #radiation
    Prad = (4./3.)*(sigma_SB/clight)*pow(T,4)
    Erad = 3*Prad/rho
    Srad = (Prad/rho + Erad)/T

    E = Egas + Erad
    P = Pgas + Prad
    S = Sgas + Srad

    beta = Pgas/P
    chiT = (Pgas + 4*Prad)/P
    chiRho = beta

    cv = 1.5*cgas/mu + 16.*(sigma_SB/clight)*pow(T,3)/rho
    cp = cv + (P*chiT*chiT)/(rho*T*chiRho)

    Gamma1 = (32-24*beta-3*beta*beta)/(24-21*beta)
    Gamma2 = (32-24*beta-3*beta*beta)/(24-18*beta - 3*beta*beta)
    Gamma3 = (32-27*beta)/(24-21*beta)

    grad_ad = (Gamma2 - 1.)/Gamma2

              #0   1     2    3    4    5    6    7       8       9     10  11   12     13     14     15
    results = [P, Prad, Pgas, E, Erad, Egas, S, Gamma1, Gamma2, Gamma3, cv, cp, beta, chiRho, chiT, grad_ad]
    return results
