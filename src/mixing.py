import numpy as np


def mix_star(star):

    # Find all convective regions
    i=0
    idx_start = []
    idx_end = []
    while i < star.N_cells:

        if star.cells_convective[i] == 1:

            idx_start.append(i)
            while i < star.N_cells and star.cells_convective[i] == 1:
                i += 1
            idx_end.append(i)

        i += 1

    N_species = len(star.X0.dtype)
    species = star.X0.dtype.names

    # Calculate instantaneous mixing
    for i in range(len(idx_start)):

        mass_species = np.zeros(len(N_species))

        # Sum up all the mass of each species in the convective region
        for j in range(N_species):
            mass_species[j] = np.sum(X0[idx_start:idx_end][species[j]])

        # Total mass for each species
        mass_total = np.sum(mass_species)

        # Normalize to get fraction of mass for each species
        for j in range(N_species):
            mass_species[j] /= mass_total

        # Set the abundances of the stars
        for k in range(idx_start, idx_end):
            for j in range(N_species):
                star.X0[k][species[j]] = mass_species[j]
