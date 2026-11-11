pars = {

    # Parameters describing the single-zone geometry
    'r_blob' :   4.62e+16, # cm
    'magfield' : 0.783, # gauss
    'lorentz' :  17.6, # Jet bulk Lorentz factor (impacts both the relativistic boosting
                      # of the external fields into the jet frame and the boost of the jet
                      # emission into the observer's frame)
    'rdiss' :    3.42e+17, # Dissipation radius in cm (cf Fig.2 of Rodrigues+2023)
                          # (impacts the density and relativistic boost of the
                          # external fields in the jet rest frame)

    # Electron injection
    # (in this case we assume a single power-law spectrum, see below)
    'elec_emin' :  2.16e1 * (const.m_e * const.c ** 2).to(u.eV).value, # eV
    'elec_emax' :  7.36e3 * (const.m_e * const.c ** 2).to(u.eV).value, # eV
    'eindex' :     2.02,
    'elum' :       1.84e+41, # erg/s

    # Proton injection
    # (in this case we assume a single power-law spectrum, see below)
    'proton_emin' :  1e2 * (const.m_p * const.c ** 2).to(u.eV).value, # eV
    'proton_emax' :  2e6 * (const.m_p * const.c ** 2).to(u.eV).value, # eV
    'pindex' :       1.0,
    'plum' :         6.29e+43, # erg/s

    # Other source parameters

    'z' : 0.19,
    'disk_lum' :  1e45 * u.erg/u.s, # accretion disk luminosity
    'black_hole_mass' : 6e8, # [solar mass units]

    'torus_temperature' : 500 * u.K,
    'blr_covering' : 0.1,
    'torus_covering' : 0.3
    }



# Below, we calculate the size of the BLR and torus and the relative Doppler
# factor of the external photons into the jet frame based on the above parameters

pars['blr_radius'] = 1e17 * (pars['disk_lum'].value / 1e45) ** 0.5  # [cm]
pars['torus_radius'] = 2.5e18 * (pars['disk_lum'].value / 1e45) ** 0.5 # [cm]
pars['blr_doppler'] = calc_doppler(pars['lorentz'], pars['blr_radius'], pars['rdiss'])
pars['torus_doppler'] = calc_doppler(pars['lorentz'],pars['torus_radius'],pars['rdiss'])

# Set homogeneous magnetic field strength
am3.set_mag_field(pars['magfield'])

# Set escape time scale (subject to escape switch - in this case escape is on)
am3.set_escape_timescale(pars['r_blob'] / const.c.cgs.value) # 4.62e16 cm == Size of the region

# Set solver time step, in this case one-tenth of the light-crossing time
am3.set_solver_time_step(1e-1 * am3.get_escape_timescale())

# Note 1:
# All the above parameters can be changed at any time during runtime. However,
# make sure to check where each parameter is used and update it accordingly.
# For example, the B-field strength is used only internally by AM3, so changing
# its value requires calling `am3.set_mag_field(new_value)`. The blob size, r_blob,
# is used internally by AM3 (for the escape rate) as well as externally (to
# calculate the particle densities in the cell below). The bulk Lorentz factor enters
# the calculation of the external fields as well as the flux in the observer's
# frame, etc.

# Note 2:
# The escape timescale would control adiabatic cooling losses as well as the reduction in particle density
# due to the expansion of the region.
# In this example, both processes are turned off because we simulate a stationary region, so there is
# no need to define it.
# am3.set_expansion_timescale(4.62e+16 / const.c.cgs.value)

# Note 3:
# The solver time step define above is quite large, and it leads to accurate results
# only for optically thin sources (e.g. most single-zone blazar models). For more
# extreme environments, it may be necessary to decrease the solver time step
# to capture the numerical hardness of the equations. As a rule of thumb, we want to
# look at the fastest process rates at the relevant particle energies.
