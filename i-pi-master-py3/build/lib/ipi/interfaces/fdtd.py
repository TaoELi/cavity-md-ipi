'''
FDTD simulation of the probe pulse on a one-dimensional grid.
'''

import numpy as np

def gaussian(z, sigma):
    return 1.0 / (np.sqrt(2.0 * np.pi) *sigma) * np.exp(-z**2 / 2.0 / sigma**2)

class fdtd1d:
    def __init__(self, dtfs=0.5, ngrids=5001, Ep_au=6e-3, omega_cminv=2320.0, rescaling_factor=1e-7):
        """Initialises one-dimensional fdtd simulation of the probe pulse.

        Args:
           dt: the default time step of FDTD simulation.
           ngrids: the number of grids of the one-dimensional spatial cell

        During simulations, a natural units are used: [c] = [hbar] = [epsilon_0] = 1
        and the time unit is fs.

        This means that I need to do unit conversion when receving data from CavMD
        au2nu and nu2au
        """
        # useful constants in units of a.u.
        self.epsilon0_au = 1.0 / 4.0 / np.pi
        self.c = 137.036
        self.mu0_au = 1.0 / self.c**2 / self.epsilon0_au
        # E-field unit conversion from natural units to atomic units and vise versa
        self.nu2au = 1.0/ 773463.5874822934
        self.au2nu = 773463.5874822934

        self.dt = dtfs
        self.dtau = dtfs * 41.341374575751 # fs2au
        self.t = 0.0
        self.ngrids = ngrids

        self.rescaling_factor = rescaling_factor

        self.nmiddle = int((ngrids-1)//2)

        # parameter for the probe pulse
        self.Ep = Ep_au * self.au2nu
        # unit converse from cm-1 to 2pi*fs-1
        self.omega_p = omega_cminv * 2.998e-5 * 2.0 * np.pi
        # unit converse from cm-1 to a.u.
        #self.omega_p = omega_cminv * 0.0000045563352812122295

        #self.dz = 2.0 * self.dtau * self.c # in units of ps-1
        self.dz = 2.0 * self.dt
        self.dtdz = 0.5 / self.c
        self.dtdzepsilon0 = 0.5 #self.dtdz / self.epsilon0_au
        self.dtdzmu0 = 0.5 #self.dtdz / self.mu0_au

        # variables for simulation
        self.Ex = np.zeros(ngrids)
        self.By = np.zeros(ngrids)

        self.cell_length = self.dz * (self.ngrids-1)
        self.Z = np.linspace(-self.cell_length / 2.0, self.cell_length / 2.0, self.ngrids)
        self.xi_x = gaussian(self.Z, sigma=self.dz)

        # parameters for absorbing boundary conditions
        self.ex_low_m1 = 0.0
        self.ex_low_m2 = 0.0
        self.ex_high_m1 = 0.0
        self.ex_high_m2 = 0.0

        self.by_low_m1 = 0.0
        self.by_low_m2 = 0.0
        self.by_high_m1 = 0.0
        self.by_high_m2 = 0.0

        # store necessary data
        self.Ex_monitor_lst = []
        self.nmonitor = int((self.ngrids - 1)//2) + 500

    def update_time(self, t0):
        self.t = t0

    def add_initial_pulse(self, fdtd_pulse_width_fs=50):
        self.sigma = fdtd_pulse_width_fs
        z_p = - self.sigma * 4
        k = self.omega_p #/ self.c
        Einitial =   self.Ep * np.sin(k * self.Z) * np.exp(-(self.Z - z_p)**2 / 2.0 / self.sigma**2)
        Binitial =  self.Ep * np.sin(k * (self.Z + self.dt * self.dz/2.)) * np.exp(-(self.Z - z_p + self.dt * self.dz/2.)**2 / 2.0 / self.sigma**2)
        self.Ex = Einitial
        self.By = Binitial

    def step(self, mux_au):
        # update the E-field
        self.Ex[1:self.ngrids] += self.dtdzepsilon0 * (self.By[0:self.ngrids-1] - self.By[1:self.ngrids])
        #Px = mux_au * self.xi_x # in atomic units
        #self.Ex -= (Px / self.epsilon0_au) * self.au2nu
        self.Ex[self.nmiddle] -= (mux_au / self.epsilon0_au) * self.au2nu * self.rescaling_factor

        # update boundary conditions
        self.Ex[0] = self.ex_low_m2
        self.ex_low_m2 = self.ex_low_m1
        self.ex_low_m1 = self.Ex[1]

        self.Ex[-1] = self.ex_high_m2
        self.ex_high_m2 = self.ex_high_m1
        self.ex_high_m1 = self.Ex[-2]

        # update the B-field
        self.By[0:self.ngrids-1] += self.dtdzmu0 * (self.Ex[0:self.ngrids-1] - self.Ex[1:self.ngrids])

        # update boundary conditions
        self.By[0] = self.by_low_m2
        self.by_low_m2 = self.by_low_m1
        self.by_low_m1 = self.By[1]

        self.By[-1] = self.by_high_m2
        self.by_high_m2 = self.by_high_m1
        self.by_high_m1 = self.By[-2]

        # finally, we also monitor the E-field value of FDTD probe pulse at a point
        Ex_monitor = self.monitor_Ex()
        # we store this value to a list
        self.Ex_monitor_lst.append(Ex_monitor)
        # oscillationlly (every 1ps by default), we output the value of the list to the folder
        if int(self.t / self.dt) % 2000 == 0:
            nsize = len(self.Ex_monitor_lst)
            data = np.zeros((nsize, 2))
            data[:,0] = np.linspace(0, self.dt*nsize-1, nsize)
            data[:,1] = np.array(self.Ex_monitor_lst)
            np.savetxt("Emonitor_fdtd.txt", data, fmt="%.4E")

    def evaluate_Ex_au(self):
        #Ex = self.dz * np.sum(self.Ex * self.xi_x)
        #return Ex * self.nu2au
        return self.Ex[self.nmiddle] * self.nu2au

    def monitor_Ex(self):
        return self.Ex[self.nmonitor] * self.nu2au

    def propagate(self, tmax=0):
        tspan = np.linspace(0.0, tmax, int(tmax/self.dt)+1)
        Ex_array = np.zeros(int(tmax/self.dt)+1)
        idx = 0
        while self.t <= tmax:
            Ex_array[idx] = self.monitor_Ex()
            self.step(mux_au=0.)
            self.update_time(self.t + self.dt)
            idx += 1
        return tspan, Ex_array

if __name__ == '__main__':
    model = fdtd1d()
    model.add_initial_pulse()
    tspan, Ex_array = model.propagate()
    # let us take fft of this data
    from scipy.fft import fft, fftfreq
    yf = fft(Ex_array)
    xf = fftfreq(np.size(tspan), 0.5) / 2.998e-5
    import matplotlib.pyplot as plt
    plt.figure()
    #plt.plot(xf, np.abs(yf))
    #plt.xlim(1500, 3000)
    #plt.plot(tspan, Ex_array)
    plt.plot(model.Ex)
    plt.show()
