import numpy as np
from scipy.integrate import odeint, solve_ivp

import constants as c
from nucleosynthesis import calc_burning
from convection import calc_star_gradT
from opacity import calc_star_opacity


class Star:

    def __init__(self, mass=1, radius=1, N_cells=100, elements=['H','He','Z']):

        # Bookkeeping
        self.N_cells = N_cells
        self.initialized = False

        # Bulk stellar properties
        self.star_mass = mass
        self.star_radius = radius
        self.star_luminosity = 0
        self.age = 0
        self.central_pressure = 0
        self.central_T = 0
        self.central_rho = 0

        # Cell properties
        self.T = np.zeros(self.N_cells)
        self.rho = np.zeros(self.N_cells)
        self.mass = np.zeros(self.N_cells)
        self.radius = np.zeros(self.N_cells)
        self.mu = np.zeros(self.N_cells)
        self.opacity = np.zeros(self.N_cells)
        self.luminosity = np.zeros(self.N_cells)

        elements_dtype = [(el,'f8') for el in elements]
        self.X0 = np.zeros(self.N_cells, dtype=elements_dtype)


        # Physics choices
        self.alpha_MLT = 2.0


    def initialize_star(self, polytropic_index):

        if polytropic_index < 0 or polytropic_index > 5-1e-5:
            print("We cannot initialize the star with this polytropic index")
            return

        def func(xi, x, n):
            theta, phi = x
            dtheta_dxi = -phi / xi**2
            dphi_dxi = theta**n * xi**2
            derivs = [dtheta_dxi, dphi_dxi]
            return derivs

        def event(xi, x, n):
            theta, phi = x
            return theta - 1.0e-10


        x0 = [1, 0]
        n = polytropic_index

        # Solve first to find maximum xi
        xi = [1e-9, 20]

        res = solve_ivp(func, xi, x0, args=(n,), events=event, atol=1.0e-20, rtol=1.0e-20)
        xi_max = res.t_events[0][0] - 1.0e-2

        # Solve again to get profile
        xi = np.linspace(1e-9, xi_max, self.N_cells)
        res = odeint(func, x0, xi, args=(n,), tfirst=True, atol=1.0e-20, rtol=1.0e-20)

        theta = res[:,0]
        phi = res[:,1]


        # Deal with any errant NaN's
        idx = np.where(np.isnan(theta))
        if np.any(idx):
            idx = idx[0][0]
            xi_max = xi[idx-1]
            theta = theta[:idx-1]
            phi = phi[:idx-1]


        print(theta)
        print(phi)




        alpha = (self.star_radius*c.Rsun) / xi_max
        central_rho = (self.star_mass*c.Msun) / (4*np.pi*alpha**3*phi[-1])
        K_const = 4*np.pi*c.G*alpha**2 / ((n+1) * central_rho**(1/n-1))

        self.radius = alpha * xi
        self.mass = 4*np.pi*alpha**3 * central_rho * phi
        self.rho = central_rho * theta**n
        self.pressure = K_const * central_rho**(1+1/n) * theta**(n+1)

        # Set the composition
        self.X0['H'] = 0.70
        self.X0['He'] = 0.28
        self.X0['Z'] = 0.02

        mu_inv = self.X0['H'] + 0.25*self.X0['He'] + 0.0625*self.X0['Z']
        self.mu = 1/mu_inv

        # Set the temperature
        self.T = self.pressure/self.rho * self.mu / c.kB * c.mp

        # Set the mass of each cell
        mass_inner = np.append([0], self.mass[:-1])
        self.dm = self.mass - mass_inner

        # Calculate the initial nuclear burning
        eps_nuc, dX_dt, dY_dt, dZ_dt = calc_burning(self.rho, self.T, self.X0)

        # Calculate the luminosity
        self.luminosity[0] = 0
        for i in range(self.N_cells-1):
            self.luminosity[i+1] = self.luminosity[i] + self.dm[i+1] * eps_nuc[i+1]


        # Set the opacity
        calc_star_opacity(self)

        # Set convective regions
        gradT, convective = calc_star_gradT(self)
        self.convective = convective


        self.initialized = True

    def evolve(self, time):

        while self.age < time:
            self.step(dt)


    def step(self, dt):

        # Calculate nuclear burning
        epsilon_total, dX_dt, dY_dt, dZ_dt = calc_burning(self.rho, self.T, self.X0)
        X0['H'] += dX_dt*dt
        X0['He'] += dY_dt*dt
        X0['Z'] += dZ_dt*dt
