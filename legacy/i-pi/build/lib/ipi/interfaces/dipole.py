"""
Given the configuration of atoms without photon modes, return the total dipole moment
"""

import numpy as np

class dipole:
    def __init__(self):
        self.have_not_set = True
        #path = os.getcwd()
        #with open(path+"/photon_params.json") as json_file:
        #    data = json.load(json_file)
        #    self.apply_self_dipole = data["add_self_dipole"]

    def set_charges(self):
        if self.have_not_set:
            self.charges = np.zeros(np.shape(self.pos[:,0]))
            # This test currently supports only H2O molecules as O H H O H H ...
            self.charges[0::3] = -0.8472
            self.charges[1::3] = 0.4236
            self.charges[2::3] = 0.4236
            print "let the charges of each atom as\n", self.charges
            self.have_not_set = False

    def update_pos(self, pos):
        self.pos = np.reshape(pos, (-1, 3))
        # Boundary condition


    def calc_dipoles_x(self):
        "calculate the x_component of the molecular dipole"
        self.set_charges()
        return self.charges * self.pos[:,0]

    def calc_dipoles_x_tot(self):
        "calculate the total value of the x_component of the molecular dipole"
        self.set_charges()
        return np.sum(self.charges * self.pos[:,0])

    def calc_dipoles_y_tot(self):
        "calculate the total value of the y_component of the molecular dipole"
        self.set_charges()
        return np.sum(self.charges * self.pos[:,1])

    def calc_dipoles_z_tot(self):
        "calculate the total value of the z_component of the molecular dipole"
        self.set_charges()
        return np.sum(self.charges * self.pos[:,2])

    def calc_dipoles_x_which(self, i=0):
        "calculate the total value of the x_component of the molecular dipole"
        self.set_charges()
        return np.sum(self.charges[i*3:i*3+3] * self.pos[i*3:i*3+3,0])

    def calc_dipoles_y_which(self, i=0):
        "calculate the total value of the y_component of the molecular dipole"
        self.set_charges()
        return np.sum(self.charges[i*3:i*3+3] * self.pos[i*3:i*3+3,1])

    def calc_dipoles_z_which(self, i=0):
        "calculate the total value of the z_component of the molecular dipole"
        self.set_charges()
        return np.sum(self.charges[i*3:i*3+3] * self.pos[i*3:i*3+3,2])

    def calc_dipoles_x_which(self, index=0):
        "calculate the total value of the x_component of the molecular dipole"
        self.set_charges()
        return np.sum(self.charges[index*3:index*3+3] * self.pos[index*3:index*3+3, 0])

    def calc_dmudx(self):
        "calculate the spatial derivative of the molecular dipole"
        return self.charges

    def calc_dmudy(self):
        "calculate the spatial derivative of the molecular dipole"
        return self.charges

    def calc_self_dipole_energy(self, pos):
        return self.charges * self.pos[:,0]

if __name__ == "__main__":
    pos = np.random.rand(18)
    print "position is\n", pos
    d = dipole()
