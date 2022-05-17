from constants import *

def eos_rhoT(rho,T,Xi):

    #composition
    X=Xi[0]
    Y=Xi[1]
    Z=Xi[2]

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

    Gamma1 = (32-24*beta-3*beta*beta)/(24-21*beta)
    Gamma2 = (32-24*beta-3*beta*beta)/(24-18*beta - 3*beta*beta)
    Gamma3 = (32-27*beta)/(24-21*beta)

    results = [P, E, S, Gamma1, Gamma2, Gamma3]
    return results
