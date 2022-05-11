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
            # Find if photon_params.json has defined an incoming pulse at initial times
            self.have_incoming_pulse = data.get("add_pulse_photon", False)
            self.have_incoming_cw = data.get("add_cw_photon", False)
            if self.have_incoming_pulse:
                self.add_pulse_direction = data.get("add_pulse_direction", 0)
                print("## Adding a pulse on photons at %d direction (0-x, 1-y, 2-z) ##" %self.add_pulse_direction)
                # pulse_params = ["E0", "tau_FWHM", "omega", "phase", "t0"]
                self.pulse_params = data.get("pulse_params", [1.0, 10.0, 3550.0, 3.14, 10.0])
                if type(self.pulse_params[2]) is list:
                    print("In fact, we add many gaussian pulses on photons at the same time...")
                    self.pulse_params[2] = np.array([x * 2.998e-5 * 2.0 * np.pi for x in self.pulse_params[2]])
                else:
                    self.pulse_params[2] *= 2.998e-5 * 2.0 * np.pi # unit converse from cm-1 to 2pi*fs-1
                    self.pulse_params[2] = np.array([self.pulse_params[2]])
                self.t = data.get("t0", 0.0)
                self.dt = data.get("dt", 0.5)
                self.transition_photon_charge = data.get("transition_photon_charge", 0.0)
                print("## add initial pulse on photons with E0 %.2E at time %.2f ##" %(self.pulse_params[0], self.pulse_params[4]))
            else:
                print("## have not set initial pulse on photons ##")
            if self.have_incoming_cw:
                self.add_cw_direction = data.get("add_cw_direction", 0)
                print("## Adding a cw field on photons at %d direction (0-x, 1-y, 2-z) ##" %self.add_cw_direction)
                # cw_params = ["E0", "omega", "phase", "tstart", "tend"]
                self.cw_params = data.get("cw_params", [1e-3, 3550.0, 3.14, 10.0, 1e4])
                # Let us add the possibility to add a few waves to the system with different frequencies
                if type(self.cw_params[0]) is list:
                    print("In fact, we add many cw waves on photons at the same time...")
                    self.have_many_cw = True
                    for idx in range(len(self.cw_params)):
                        self.cw_params[idx][1] *= 2.998e-5 * 2.0 * np.pi
                        print("## For cw No.%d, add initial cw with E0 %.2E at time %.2f ending at \
                            %.2f##" %(idx, self.cw_params[idx][0], self.cw_params[idx][3], self.cw_params[idx][4]) )
                else:
                    self.have_many_cw = False
                    self.cw_params[1] *= 2.998e-5 * 2.0 * np.pi # unit converse from cm-1 to 2pi*fs-1
                    print("## add initial cw with E0 %.2E at time %.2f ending at \
                        %.2f##" %(self.cw_params[0], self.cw_params[3], self.cw_params[4]) )
                self.t = data.get("t0", 0.0)
                self.dt = data.get("dt", 0.5)
                self.transition_photon_charge = data.get("transition_photon_charge", 0.0)
            else:
                print("## have not set initial cw on photons ##")
        else:
            print("No photon applied")
        print("### End of initialization ###")

    def update_pos(self, pos):
        #pos = [x1, y1, z1, x2, y2, z2, ...; x1, y1, z1, x2, y2, z2, ...]
        self.pos = pos

    '''
    def add_pulse(self, mf):
        if self.have_incoming_pulse:
            self.t += self.dt
            if self.t > self.pulse_params[4] and self.t < self.pulse_params[4] + self.pulse_params[1]*8.0:
                t = self.t - self.pulse_params[4] - self.pulse_params[1]*4.0
                Ex = np.sum(self.pulse_params[0] * np.exp(-t**2 / self.pulse_params[1]**2 \
                    * 2.0 * np.log(2.0)) * np.cos(self.pulse_params[2]*self.t + self.pulse_params[3]))
                if self.add_pulse_direction == 0:
                    mf[0:self.nmodes*3:3] += Ex * self.transition_photon_charge
                elif self.add_pulse_direction == 1:
                    mf[self.nmodes*3+1::3] += Ex * self.transition_photon_charge
    '''
    def add_pulse(self, mf):
        if self.have_incoming_pulse:
            self.t += self.dt
            if self.t > self.pulse_params[4] and self.t < self.pulse_params[4] + self.pulse_params[1]*4.0:
                t = self.t - self.pulse_params[4] #- self.pulse_params[1]*4.0
                Ex = np.sum(self.pulse_params[0] * np.exp(-t**2 / self.pulse_params[1]**2 \
                    * 2.0 * np.log(2.0)) * np.sin(self.pulse_params[2]*self.t + self.pulse_params[3]))
                if self.add_pulse_direction == 0:
                    mf[0:self.nmodes*3:3] += Ex * self.transition_photon_charge
                elif self.add_pulse_direction == 1:
                    mf[self.nmodes*3+1::3] += Ex * self.transition_photon_charge

    def add_cw(self, mf, phase):
        if self.have_incoming_cw:
            self.t += self.dt
            if not self.have_many_cw:
                # cw is not initalized for the molecules and we use random phase
                if phase is None and self.cw_params[2] == "RANDOM_PHASE":
                    # we reset the pulse phase
                    self.cw_params[2] = np.random.rand() * 2.0 * np.pi
                    print("Initialize the cw phase for photons as", self.cw_params[2])
                # cw is not initalized for the molecules and we do not use random phase; do nothing
                # cw is initalized for the molecules and we use random phase; assign phase to our parameters
                # cw is initalized for the molecules and we do not use random phase; do noting
                if phase is not None and self.cw_params[2] == "RANDOM_PHASE":
                    self.cw_params[2] = phase

                if self.t > self.cw_params[3] and self.t < self.cw_params[4]:
                    Ex = self.cw_params[0] * np.cos(self.cw_params[1]*self.t + self.cw_params[2])
                    if self.add_cw_direction == 0:
                        mf[0:self.nmodes*3:3] += Ex * self.transition_photon_charge
                    elif self.add_cw_direction == 1:
                        mf[self.nmodes*3+1::3] += Ex * self.transition_photon_charge
            else:
                if phase is None:
                    for cw_params in self.cw_params:
                        if cw_params[2] == "RANDOM_PHASE":
                            cw_params[2] = np.random.rand() * 2.0 * np.pi
                            print("Initialize the cw phase for photons as", cw_params[2])
                else:
                    for cw_params, mphase in zip(self.cw_params, phase):
                        if cw_params[2] == "RANDOM_PHASE":
                            cw_params[2] = mphase
                Ex = 0.0
                for cw_params in self.cw_params:
                    if self.t > cw_params[3] and self.t < cw_params[4]:
                        Ex += np.sum(cw_params[0] * np.cos(cw_params[1]*self.t + cw_params[2]))
                if self.add_cw_direction == 0:
                    mf[0:self.nmodes*3:3] += Ex * self.transition_photon_charge
                elif self.add_cw_direction == 1:
                    mf[self.nmodes*3+1::3] += Ex * self.transition_photon_charge


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
