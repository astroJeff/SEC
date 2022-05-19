import numpy as np

def calc_burning(rho, T, X0):

    # Energy production rate in erg/s/g
    epsilon_pp = calc_pp(rho, T, X0)
    epsilon_cno = calc_cno(rho, T, X0)
    epsilon_triple_alpha = calc_triple_alpha(rho, T, X0)
    epsilon_total = epsilon_pp + epsilon_cno + epsilon_triple_alpha

    # Energy released per reaction in erg/g
    E_pp = 6.3e18
    E_cno = 6.0e18
    E_triple_alpha = 6.0e17

    # Change compositions
    dX_dt = -epsilon_pp/E_pp - epsilon_cno/E_cno
    dY_dt = epsilon_pp/E_pp + epsilon_cno/E_cno - epsilon_triple_alpha/E_triple_alpha
    dZ_dt = epsilon_triple_alpha/E_triple_alpha

    return epsilon_total, dX_dt, dY_dt, dZ_dt


def calc_pp(rho, T, X0):

    X = X0['H']
    Y = X0['He']
    Z = X0['Z']

    T6 = T/1e6

    f1 = 1 + 0.25*np.sqrt(rho) * T6**(-3/2)
    g1 = 1 + 0.0012*T6**(1/3) + 0.0078*T6**(2/3) + 0.0006*T6

    epsilon_pp = 2.06e6 * f1 * g1 * X**2*rho * T6**(-2/3) * np.exp(-33.81 * T6**(-1/3))

    return epsilon_pp


def calc_cno(rho, T, X0):

    X = X0['H']
    Y = X0['He']
    Z = X0['Z']

    T6 = T/1e6

    f1 = 1 + 1.75*np.sqrt(rho) * T6**(-3/2)
    g1 = 1 + 0.0027*T6**(1/3) - 0.0037*T6**(2/3) - 0.0007*T6

    epsilon_cno = 7.94e27 * f1 * g1 * rho * X*Z * T6**(-2/3) * np.exp(-152.313 * T6**(-1/3))

    return epsilon_cno

def calc_triple_alpha(rho, T, X0):

    X = X0['H']
    Y = X0['He']
    Z = X0['Z']

    T8 = T/1e8

    f1 = np.exp(2.4e-3 * np.sqrt(rho) * T8**(-3/2))

    epsilon_triple_alpha = 3.46e11 * f1 * Y**3 * rho**2 * T8**(-3) * np.exp(-43.2 / T8)

    return epsilon_triple_alpha
