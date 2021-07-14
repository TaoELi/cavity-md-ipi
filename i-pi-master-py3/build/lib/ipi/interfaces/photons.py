"""
Contain the information for photon modes
"""

import numpy as np
import os
import json

class photon:
    def __init__(self):
        self.hartree_to_cm_minus1 = 219474.63
        print("### Initialise a photon mode from photon_params.json ###")
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
            print("set photon with effective mass %.3f a.u." %self.mass)
            print("set photon freq as %.3f cm-1" %(self.hartree_to_cm_minus1 * self.freq))
            print("set photon E-field magnitude as %.3E" %self.E0)
            self.pos = np.zeros(3, np.float64)
            self.coeff_self = self.E0**2 / self.mass / self.freq**2
        else:
            print("No photon applied")
        print("### End of initialization ###")

    def update_pos(self, pos):
        self.pos = pos.copy()

    def obtain_potential(self):
        return self.pot_coeff * np.sum(self.pos**2)

    def obtain_Ex(self):
        return self.E0 * self.pos[0]

    def obtain_Ey(self):
        return self.E0 * self.pos[1]

    def obtain_dEdx(self):
        return self.E0

    def obtain_bare_force(self):
        return - (2.0 * self.pot_coeff) * self.pos

class photons:
    def __init__(self):
        self.hartree_to_cm_minus1 = 219474.63
        path = os.getcwd()
        self.apply_photon = True
        try:
            with open(path+"/photon_params.json") as json_file:
                data = json.load(json_file)
        except:
            print("Not found photon_params.json, do conventional nuclear dynamics")
            self.apply_photon = False
        if self.apply_photon:
            self.apply_photon = data["apply_photon"]
        if self.apply_photon:
            print("### Initialise many photon modes from photon_params.json ###")
            self.nmodes = data.get("n_modes", 1)
            self.nphoton = self.nmodes * 2
            self.mass = data["eff_mass"]
            # Representing the freq and E0 for the fundamental two modes (in x and y directions)
            self.freq = data["freqs_cm"] / self.hartree_to_cm_minus1
            self.E0 = data["E0"]
            print("set fundamental photon with effective mass %.3f a.u." %self.mass)
            print("set fundamental photon freq as %.3f cm-1" %(self.hartree_to_cm_minus1 * self.freq))
            print("set fundamental photon E-field magnitude as %.3E" %self.E0)
            print("The system has %d * 2 photon modes" %(self.nmodes))
            # Obtaining the freq and E0 list for N*2 modes (in x and y directions) [x1, x2, x3, ..., xN, y1, y2, ...]
            self.freq_lst = self.freq * np.array([n for n in range(1, self.nmodes+1)]*2)
            self.E0_lst = self.E0 * np.ones(self.nmodes*2) * self.freq_lst / self.freq
            print("set photon freqs as [cm-1]:")
            print(self.hartree_to_cm_minus1 * self.freq_lst)
            print("set photon E-field magnitudes as")
            print(self.E0_lst)
            self.pos = np.zeros(self.nphoton*3, np.float64)
            self.coeff_self = np.sum(self.E0_lst[0:self.nmodes]**2 / self.mass / self.freq_lst[0:self.nmodes]**2)
            self.pot_coeff = 0.5 * self.mass * self.freq_lst**2
            self.pot_coeff2 = self.mass * np.reshape(np.array([[x,x,x] for x in self.freq_lst]), -1)**2
        else:
            print("No photon applied")
        print("### End of initialization ###")

    def update_pos(self, pos):
        #pos = [x1, y1, z1, x2, y2, z2, ...; x1, y1, z1, x2, y2, z2, ...]
        self.pos = pos

    def obtain_potential(self):
        return np.sum(self.pot_coeff2 * self.pos**2) / 2.0

    def obtain_Ex(self):
        return np.sum(self.E0_lst[0:self.nmodes] * self.pos[0:self.nmodes*3:3])

    def obtain_Ey(self):
        return np.sum(self.E0_lst[self.nmodes:] * self.pos[self.nmodes*3+1::3])

    def calc_photon_force(self, mu_x, mu_y):
        f = - self.pot_coeff2 * self.pos
        f[0:self.nmodes*3:3] -= mu_x * self.E0_lst[0:self.nmodes]
        f[self.nmodes*3+1::3] -=  mu_y * self.E0_lst[self.nmodes:]
        return f


if __name__ == "__main__":
    p = photons()
