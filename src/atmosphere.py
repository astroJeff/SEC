import constants as c


def calc_outer_BC(star):

    g_grav = -c.G * star.star_mass / star.star_radius**2

    # For opacity
    zion = [1,2,8]
    aion = [1,4,16]
    ionmax = 3
    pep_in, xne_in, eta_in = 0, 0, 0

    kappa = calc_opacity(star.T[-1],
                         star.rho[-1],
                         star.X0[-1],
                         zion,
                         aion,
                         ionmax,
                         pep_in,
                         xne_in,
                         eta_in)

    # Simple photosphere
    pressure_outer = 2/3 * g_grav / kappa
    T_outer = (star.star_luminosity / (4*np.pi*c.sigma_SB*star.star_radius**2))**(1/4)

    return pressure_outer, T_outer
