"""
Contain the information for photon modes
"""

import numpy as np
import os
import json

class photon:
    def __init__(self):
        self.hartree_to_cm_minus1 = 220000.0
        print "### Initialise a photon mode from photon_params.json ###"
        path = os.getcwd()
        with open(path+"/photon_params.json") as json_file:
            data = json.load(json_file)
        self.apply_photon = data["apply_photon"]
        if self.apply_photon:
            self.mass = data["eff_mass"]
            self.freq = data["freqs_cm"] / self.hartree_to_cm_minus1
            self.pot_coeff = 0.5 * self.mass * self.freq**2
            self.E0 = data["E0"]
            self.pos = np.zeros(3, np.float64)
            print "set photon with effective mass %.3f a.u." %self.mass
            print "set photon freq as %.3f cm-1" %(self.hartree_to_cm_minus1 * self.freq)
            print "set photon E-field magnitude as %.3E" %self.E0
            self.pos = np.zeros(3, np.float64)
        else:
            print "No photon applied"
        print "### End of initialization ###"

    def update_pos(self, pos):
        self.pos = pos.copy()
        #print "position is", self.pos

    def obtain_potential(self):
        return self.pot_coeff * np.sum(self.pos**2)

    def obtain_Ex(self):
        return self.E0 * self.pos[0]

    def obtain_dEdx(self):
        return self.E0

    def obtain_bare_force(self):
        return - 2.0 * self.pot_coeff * self.pos

if __name__ == "__main__":
    pos = np.random.rand(18)
    print "position is\n", pos
    p = photon()
    print p.apply_photon
    p.update_pos(pos[-3:])
    print p.obtain_potential()
    print p.obtain_Ex()
    print p.obtain_dEdx()
    print p.obtain_bare_force()
