# Create electron spectrum, in this case a simple power law.
# The same can also be achieved using the template powerlaw function:
# am3.set_electron_powerlaw_injection_parameters_parameters((elum, elec_emin, elec_emin, elec_emax, eindex, eindex, 1.0))

egrid = am3.get_egrid_lep()
epowerlaw = (egrid ** (2 - pars['eindex'])
                * (egrid >= pars['elec_emin'])
                * np.exp(- egrid / pars['elec_emax'])
            )
# Integrate spectrum in erg/s
etrapz = np.trapezoid(epowerlaw / egrid, egrid)

# Normalize it to the desired electron lumiosity
enormalized = epowerlaw * pars['elum'] / etrapz # erg/s

# Convert it to an energy density injection rate
volume = 4/3 * np.pi * pars['r_blob'] ** 3
enormalized /= volume # erg/cm3/s

# And finally to a particle density injection rate
enormalized /= egrid * u.eV.to(u.erg) # cm-3.s-1



# Create proton spectrum, in this case also a simple power law.
# The same can be achieved using the template powerlaw function:
# io.set_proton_powerlaw_injection_parameters_parameters((plum, proton_emin, proton_emin, proton_emax, pindex, pindex, 1.0))

pgrid = am3.get_egrid_had()
ppowerlaw = (pgrid ** (2 - pars['pindex'])
                * (pgrid >= pars['proton_emin'])
                * np.exp(- pgrid / pars['proton_emax'])
            )

# Integrate spectrum in erg/s
ptrapz = np.trapezoid(ppowerlaw / pgrid, pgrid)

# Normalize it to the desired proton lumiosity
pnormalized = ppowerlaw * pars['plum'] / ptrapz # erg/s

# Convert it to an energy density injection rate
volume = 4/3 * np.pi * pars['r_blob'] ** 3
pnormalized /= volume # erg/cm3/s

# And finally to a particle density injection rate
pnormalized /= pgrid * u.eV.to(u.erg) # cm-3.s-1
