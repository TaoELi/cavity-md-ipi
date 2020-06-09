"""
Given the configuration of atoms without photon modes, return the total dipole moment
"""

import numpy as np

class dipole:
    def __init__(self):
        self.have_not_set = True

    def set_charges(self, pos):
        if self.have_not_set:
            self.charges = np.zeros(int(len(pos) / 3), np.float64)
            # This test currently supports only H2O molecules as O H H O H H ...
            self.charges[0::3] = -1.04
            self.charges[1::3] = 0.52
            self.charges[2::3] = 0.52
            print "let the charges of each atom as\n", self.charges
            self.have_not_set = False

    def calc_dipoles_x(self, pos):
        "calculate the x_component of the molecular dipole"
        self.set_charges(pos)
        return self.charges * pos[0::3]

    def calc_dmudx(self, pos):
        "calculate the spatial derivative of the molecular dipole"
        return self.charges

if __name__ == "__main__":
    pos = np.random.rand(18)
    print "position is\n", pos
    d = dipole()
    print d.calc_dipoles_x(pos)
    print d.calc_dipoles_x(pos)
    print d.calc_dmudx(pos)
