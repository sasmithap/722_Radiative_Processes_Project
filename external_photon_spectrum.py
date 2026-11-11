# AM3 photon grid
egrid_jetframe = am3.get_egrid_photons() * u.eV

# Set up array for adding up external fields
external_photons = np.zeros(egrid_jetframe.size) * u.GeV / u.cm**3

# Scattered thermal disk emission
disky = ShakuraFlux(egrid_jetframe / pars['blr_doppler'],
                    pars['disk_lum'],
                    pars['black_hole_mass'],
                    3.0) # erg/s, black hole frame

# Broad line emission
hybl = BroadLine(egrid_jetframe / pars['blr_doppler'],
                    10.2*u.eV, 10.2*u.eV/20, # H Ly alpha line energy and width
                    pars['disk_lum'] * pars['blr_covering']) # [erg/s]
hebl = BroadLine(egrid_jetframe / pars['blr_doppler'],
                    40.8*u.eV, 40.8*u.eV/20, # He Ly alpha line energy and width
                    pars['disk_lum'] * pars['blr_covering'] * 0.5) # [erg/s]
# Convert BLR components to jet frame
blr_to_jet = convert_lum_to_density_in_jet(pars['rdiss'],
                                           pars['lorentz'],
                                           pars['blr_radius']) # erg/s -> GeV/cm3
blr_to_jet *= u.GeV / u.cm ** 3 / (u.erg / u.s) # give it units

blr_pho = (disky * 0.01 + hybl + hebl) * blr_to_jet # GeV/cm3

external_photons += blr_pho

# Dust torus
torusy = PlanckDistribution(egrid_jetframe / pars['torus_doppler'],
                            pars['torus_temperature'],
                            pars['disk_lum'] * pars['torus_covering']) # [erg/s]
# Convert torus emission to jet frame
tor_to_jet = convert_lum_to_density_in_jet(pars['rdiss'],
                                           pars['lorentz'],
                                           pars['torus_radius']) # erg/s -> GeV/cm3
tor_to_jet *= u.GeV / u.cm**3 / (u.erg/u.s)  # give it units
# Add torus to BLR components
torus_pho = torusy * tor_to_jet
external_photons += torus_pho

# Convert summed up components from energy density to photon density in jet frame
external_photonspectrum = (external_photons / egrid_jetframe).to(u.cm ** -3).value # cm-3

# Convert photon density to density ijnjection rate
external_photonspectrum /= am3.get_escape_timescale() # cm-3.s-1
