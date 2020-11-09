import json
import numpy as np

single_charges = [0.6512, -0.3256, -0.3256]

data = {}
data["apply_photon"] = True
#data["n_modes"] = 1
data["eff_mass"] = 1.0
data["freqs_cm"] = 2320.0
data["E0"] = 0.0005

# cw parameters
data["add_cw"] = True
data["add_cw_direction"] = 0
data["cw_params"] = [0.05, 2241.0, "RANDOM_PHASE", 100.0, 600.0]
data["cw_atoms"] = [-1]
data["dt"] = 0.5

charges = single_charges*343




print("charges length = ", len(charges))
# I need partial charge at each atom (C, O, O)
charges = [float("%.5f" %x) for x in charges]

data["charge_array"] = charges

with open('test.json', 'w') as outfile:
    json.dump(data, outfile)

