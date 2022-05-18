from constants import *
from eos import eos_rhoT
from opacity import calc_opacity
from numpy import sqrt

#Henyey MLT
mlt_f1 = 0.125
mlt_f2 = 0.5
mlt_f3 = 24.
mlt_f4 = 3.

def gradT_MLT(rho,T,Xi,grad_rad,grav,alpha_MLT):
    eos=eos_rhoT(rho,T,Xi)
    P = eos[0]
    cp = eos[11]
    Q0 = eos[14]/eos[13] #chiT/chiRho
    grad_ad = eos[15]

    opac=calc_opacity(T,rho,Xi,[1,2,8],[1,4,16],3,0.0, 0.0, 0.0)
    kappa=opac[2]

    
    pressure_scale_height = P/(grav*rho)
    Lambda = alpha_MLT * pressure_scale_height
    
    omega = Lambda * rho * kappa
    a0 = 0.1875 * mlt_f2 * mlt_f3/(1. + mlt_f4/(omega*omega))
    A_0 = sqrt(mlt_f1*P*Q0*rho)
    A_1 = 4.* A_0 * cp
    A_2 = alpha_MLT*omega*(1. + mlt_f4/(omega*omega))
    A = (A_1*A_2)/(mlt_f3*crad*clight*pow(T,3))

    B3 = (A*A/a0)*(grad_rad-grad_ad)
    f = -2 + 9*a0 + 27*pow(a0,3)*B3

    if f > 1e100:
        f0=f
    else:
        f0 = f*f + 4*pow(-1 + 3*a0,3)
        if f0 < 0:
            f0=f
        else:
            f0=sqrt(f0)

    f1 = -2 + 9*a0 + 27*pow(a0,3)*B3 + f0
    f1 = pow(f1,1./3.)
    f2 = 2*pow(2,1./3.)*(1 - 3*a0)/f1

    Gamma = (pow(4,1./3.)*f1 + f2 -2)/(6*a0)
    Zeta = min(1.,max(0.,pow(Gamma,3)/B3))

    gradT = (1 - Zeta)*grad_rad + Zeta*grad_ad

    return gradT
