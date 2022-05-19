def calc_radius(star):

    # Set the radius
    ln_radius = np.zeros(N_cells)
    for i in range(star.N_cells-1):
        ln_radius[i] = 1/3 * np.log(r[i+1]**3 + 3/4/np.pi * star.dm/star.rho)


    # Set the pressure
    for k in range(star.N_cells-1):
        m_bar = (star.dm[k] + star.dm[k-1]) / 2

        P[k-1] = P[k] + m_bar * (-c.G*star.mass[k]/(4*np.pi*star.radius[k]**4))

    # Set the temperature
    gradT = calc_star_gradT(star)
    for T in range(star.N_cells-1):
        m_bar = (star.dm[k] + star.dm[k-1]) / 2
        T_bar = (star.T[k] + star.T[k-1]) / 2
        P_bar = (star.pressure[k] + star.pressure[k-1]) / 2

        T[k-1] = T[k] + m_bar * (gradT * -c.G*star.mass[k]/(4*np.pi*star.radius[k]**4) * T_bar/P_bar)

    # Energy equation
    for k in range(star.N_cells-1):
        L[k] = L[k+1] + star.dm * (eps_nuc + eps_grav)






    pressure_at_surface, temperature_at_surface = calc_outer_BC(star)
