import numpy as np
from scipy.integrate import odeint, solve_ivp

import constants as c
from nucleosynthesis import calc_burning

class Star:

    def __init__(self, mass=1, radius=1, N_cells=100, elements=['H','He','Z']):

        # Bookkeeping
        self.N_cells = N_cells
        self.initialized = False

        # Bulk stellar properties
        self.mass = mass
        self.radius = radius
        self.luminosity = 0
        self.age = 0
        self.central_pressure = 0
        self.central_T = 0
        self.central_rho = 0

        # Cell properties
        self.cells_T = np.zeros(self.N_cells)
        self.cells_rho = np.zeros(self.N_cells)
        self.cells_mass = np.zeros(self.N_cells)
        self.cells_radius = np.zeros(self.N_cells)
        self.cells_mu = np.zeros(self.N_cells)

        elements_dtype = [(el,'f8') for el in elements]
        self.X0 = np.zeros(self.N_cells, dtype=elements_dtype)


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
        res = odeint(func, x0, xi, args=(n,), tfirst=True)

        theta = res[:,0]
        phi = res[:,1]

        alpha = (self.radius*c.Rsun) / xi_max
        central_rho = (self.mass*c.Msun) / (4*np.pi*alpha**3*phi[-1])
        K_const = 4*np.pi*c.G*alpha**2 / ((n+1) * central_rho**(1/n-1))

        self.cells_radius = alpha * xi
        self.cells_mass = 4*np.pi*alpha**3 * central_rho * phi
        self.cells_density = central_rho * theta**n
        self.cells_pressure = K_const * central_rho**(1+1/n) * theta**(n+1)

        # Set the composition
        self.X0['H'] = 0.70
        self.X0['He'] = 0.28
        self.X0['Z'] = 0.02

        mu_inv = self.X0['H'] + 0.25*self.X0['He'] + 0.0625*self.X0['Z']
        self.cells_mu = 1/mu_inv

        # Set the temperature
        self.cells_T = self.cells_pressure/self.cells_density * self.cells_mu / c.kB * c.mp

        self.initialized = True

    def evolve(self, time):

        while self.age < time:
            self.step()


    def step(self, dt):

        # Calculate nuclear burning
        epsilon_total, dX_dt, dY_dt, dZ_dt = calc_burning(self.rho, self.T, self.X0)
        X0['H'] += dX_dt*dt
        X0['He'] += dY_dt*dt
        X0['Z'] += dZ_dt*dt
