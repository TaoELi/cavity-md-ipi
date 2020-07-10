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
            self.apply_inplane_dist = data.get("apply_inplane_dist", False)
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
            if not self.apply_inplane_dist:
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
                # We sample many inplane wave vectors in 2d geometry, resulting in a polariton bath!
                # First, we read controlling paramters to define the inplane wave vectors
                # note that the total number of cavity modes is still self.nmodes
                self.kparallel_max = data.get("inplane_kmax_ratio", 0.2) * self.freq
                # The frequency of cavity mode is calculated by sqrt(freq**2 + kparallel**2)
                # We need to define discrete k points in a square
                kxky_array = np.zeros((2, self.nmodes))
                npoints_1d = int(np.sqrt(self.nmodes))
                if npoints_1d is 1:
                    dk = 0.0
                else:
                    dk = self.kparallel_max * 2.0 / (npoints_1d - 1)
                num = 0
                for i in range(npoints_1d):
                    for j in range(npoints_1d):
                        kxky_array[0, num] = dk * i - self.kparallel_max
                        kxky_array[1, num] = dk * j - self.kparallel_max
                        num += 1
                if self.nmodes is 1:
                    kxky_array[0,0] = 0.0
                    kxky_array[1,0] = 0.0
                print("# Inplane wave vector lattice is (kx, ky) [in cm-1] #")
                print(kxky_array * self.hartree_to_cm_minus1)
                # Now, we calculate frequency list
                self.freq_lst = np.zeros(self.nmodes*2)
                for i in range(self.nmodes):
                    self.freq_lst[i] = self.freq_lst[i+self.nmodes] = np.sqrt(np.sum(kxky_array[:, i]**2) + self.freq**2)
                print("# Photon frequencies are [in cm-1] #")
                print(self.freq_lst * self.hartree_to_cm_minus1)
                # Calculate the vxi unit vectors for each cavity mode
                self.vxi_1 = np.zeros((self.nmodes, 3))
                self.vxi_2 = np.zeros((self.nmodes, 3))
                self.vxi_tot = np.zeros((self.nmodes*2, 3))
                for i in range(self.nmodes):
                    beta = np.arccos(self.freq / self.freq_lst[i])
                    alpha = np.arctan2(kxky_array[1, i], kxky_array[0, i])
                    self.vxi_1[i,:] = np.array([-np.sin(alpha), np.cos(alpha), 0.0])
                    self.vxi_2[i,:] = np.array([np.cos(alpha) * np.cos(beta), np.sin(alpha) * np.cos(beta), -np.sin(beta)])
                self.vxi_tot[0:self.nmodes,:] = self.vxi_1
                self.vxi_tot[self.nmodes:,:] = self.vxi_2
                print("# Unit vectors of photon polarization are #")
                print(self.vxi_1)
                print(self.vxi_2)
                print("#combining together#")
                print(self.vxi_tot)
                # Calculate the E0 lst
                self.E0_lst = self.E0 * np.ones(self.nmodes*2) * self.freq_lst / self.freq
                print("set photon E-field magnitudes as")
                print(self.E0_lst)
                self.pos = np.zeros(self.nphoton*3, np.float64)
                self.coeff_self_array = self.E0_lst**2 / self.mass / self.freq_lst**2
                self.pot_coeff = 0.5 * self.mass * self.freq_lst**2
                self.pot_coeff2 = self.mass * np.reshape(np.array([[x,x,x] for x in self.freq_lst]), -1)**2
                print("pot_coeff2 is")
                print(np.sqrt(self.pot_coeff2.reshape(-1,3) / self.mass * self.hartree_to_cm_minus1**2))
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

    def calc_dipole_projecting_vxi(self, dipole_tot_array):
        ''' dipole_tot_array dimension 3x1, return nphoton different dtot,klambda '''
        return np.sum(dipole_tot_array * self.vxi_tot, axis=-1)

    def calc_photon_force(self, mu_x, mu_y):
        f = - self.pot_coeff2 * self.pos
        f[0:self.nmodes*3:3] -= mu_x * self.E0_lst[0:self.nmodes]
        f[self.nmodes*3+1::3] -=  mu_y * self.E0_lst[self.nmodes:]
        return f

    def calc_photon_force_inplane(self, dipole_tot_array):
        f = - self.pot_coeff2 * self.pos
        dtot_klambda = self.calc_dipole_projecting_vxi(dipole_tot_array)
        f[0::3] -= dtot_klambda * self.E0_lst # This seems to be problematic
        #print np.reshape(f,(-1,3))
        return f

    def calc_nuclear_force_cavity(self, dipole_tot_array, charge_array):
        dtot_klambda = self.calc_dipole_projecting_vxi(dipole_tot_array)
        #print "dtot_klambda is"
        #print dtot_klambda
        inbracket_klambda_array = -self.E0_lst * self.pos[0::3] - self.coeff_self_array * dtot_klambda
        #print "inbracket_klambda_array is"
        #print inbracket_klambda_array
        photon_coeff = np.sum(np.reshape(inbracket_klambda_array, (-1,1)) * self.vxi_tot, axis=0)
        #print "photon_coeff before reshape is"
        #print photon_coeff
        photon_coeff = photon_coeff.reshape(1, 3)
        #print "photon_coeff after reshape is"
        #print photon_coeff
        charge_array = np.reshape(charge_array, (-1, 1))
        #print "charge array after reshape is"
        #print charge_array
        result = np.kron(charge_array, photon_coeff)
        #print "product between charge array and photon coeff is"
        #print result
        #print "Finally, we need", np.shape(result.reshape(-1))
        return result.reshape(-1)


if __name__ == "__main__":
    p = photons()
    x = np.reshape(np.array([0.0, 1.0, 0.0]), (1,3))
    y = np.array([0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    #print p.calc_dipole_projecting_vxi(x)
    #print p.calc_photon_force_inplane(x)
    print(np.reshape(p.calc_nuclear_force_cavity(x, y), (-1,3)))
