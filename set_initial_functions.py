%%time

# Initialize the AM3 module
am3 = am3.AM3()

# Adjust min photon energy, min neutrino energy, and max energy for all particles [eV].
# This allows us to optimize the performance depending on the needs of the specific science case
am3.update_energy_grid(1e-6, 1e3, 1e18)

# Manually define max proton and electron energy
am3.set_estimate_max_energies(0)

# Turn on for keeping track of the different SED components, for plotting purposes;
# turn off for efficiency.
am3.set_process_parse_sed(1)

# Hadronic processes on
am3.set_process_hadronic(1)

# Keep track of positrons and electrons separately
am3.set_process_merge_positrons_into_electrons(0)

# Escape (see documentation)
am3.set_process_escape(1)

# Expansion of radiation zone (see documentation). This one is stationary.
am3.set_process_expansion(0)
am3.set_process_adiabatic_cooling(0)

# Electron synchrotron and synchrotron self-absorption
am3.set_process_electron_syn(1)
am3.set_process_ssa(1)

# Proton synchrotron - subdominant in this case, but may be relevant for B>~10 G
am3.set_process_proton_syn(1)

# Quantum synchrotron for pairs off by default.
# In AGN simulations this effect can typically be neglected, so turn off for efficiency
am3.set_process_quantum_syn(0)

# Inverse Compton by electrons and protons
am3.set_process_electron_compton(1)
am3.set_process_proton_compton(1)

# Direct Compton turned off by default
am3.set_process_compton_photon_energy_loss(0)


# Synchrotron and inverse Compton by muons and pions - off by default.
# In AGN simulations this effect can typically be neglected, so turn off for efficiency
am3.set_process_muon_syn(0)
am3.set_process_pion_syn(0)
am3.set_process_muon_compton(0)
am3.set_process_pion_compton(0)

# Secondary particle decay
am3.set_process_pion_decay(1)
am3.set_process_muon_decay(1)

# Photon annihilation (gamma gamma -> e- e+)
am3.set_process_annihilation(1)
am3.set_optimize_annihilation_pair_emission(1) # optimization

# Bethe-Heitler pair production (p gamma -> p e+ e-)
am3.set_process_bethe_heitler(1)
am3.set_optimize_bethe_heitler_outgoing_pairs_grid(1)
am3.set_optimize_bethe_heitler_incoming_protons_min(1e12)
am3.set_optimize_bethe_heitler_target_photon_max(1e6)

# Photo-pion production (nucleon gamma -> nucleon pion)
am3.set_process_photopion(1)
am3.set_optimize_photopion_target_photon_grid(1)
am3.set_optimize_photopion_target_photon_max(1e6)

am3.set_profile_timing(1) # Used for profiling, as in the bottom plot. Turn off for efficiency.

# Initialize the AM3 kernels with the above settings
am3.init_kernels()
