# This script is used to capture many important data from each xyz trajectory

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import os, sys, time
import math
from itertools import islice
from scipy import fftpack
import glob
import json

def smooth(x,window_len=11,window='hamming'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """
    if window_len<3:
        return x

    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[window_len//2-1:-window_len//2]


class MD_Analysis:
    def __init__(self, xyz_filename, dtfs = 2, nframe_max=50001):
        self.xyz_filename = xyz_filename
        self.dtfs = dtfs
        self.dtau = dtfs * 1e-15 / 2.418884326e-17
        self.nframe_max = nframe_max
        self.natoms = 0
        self.labels = []
        self.traj = []
        self.load_xyz(self.xyz_filename)
        # After read xyz file, now we calculate different properties

    def read_xyz_info(self, xyz_filename):
        print("Reading info of atoms")
        next_n_lines = lambda file_opened, N : [x.strip() for x in islice(file_opened, N)]
        self.labels = []
        with open(xyz_filename, 'r') as myfile:
            self.natoms = int(myfile.readline().strip())
            myfile.readline()
            data = next_n_lines(myfile, self.natoms)
            for i in range(self.natoms):
                self.labels.append( data[i].split()[0] )

    def load_xyz(self, xyz_filename):
        # Load info of how many water molecules are available
        self.read_xyz_info(xyz_filename)
        # Load trajectory
        try:
            print("Try to load file from %s.npy" %xyz_filename)
            self.traj = np.load("%s.npy" %xyz_filename)
            print("Loaded file from %s.npy" %xyz_filename)
        except:
            next_n_data = lambda file_opened, N : [x.strip().split()[1:] for x in islice(file_opened, N)]
            self.traj = np.zeros((self.natoms, 3, self.nframe_max))
            print("Reading trajs ...")
            with open(xyz_filename, 'r') as myfile:
                for n in range(self.nframe_max):
                    try:
                        data = next_n_data(myfile, self.natoms+2)[2:]
                        coords = np.asarray(data).astype(np.float)
                        self.traj[:, :, n] = coords
                        if (n % 1000 == 0):
                            print("Finish No.%d frame" %n)
                    except:
                        print("Read No.%d frame but has error" %n)
                        break
            # Save data
            np.save(xyz_filename + ".npy", self.traj)

    def cacl_center_of_mass_water_traj(self):
        print("Calculating Center of Mass traj of water")
        self.nwaters = self.natoms // 3
        mO, mH, mD = 15.9994, 1.00794, 2.0141
        O_traj = self.traj[0:self.nwaters*3:3, :, :]
        H1_traj = self.traj[1:self.nwaters*3:3, :, :]
        H2_traj = self.traj[2:self.nwaters*3:3, :, :]
        self.com_water_traj = (O_traj * mO + H1_traj * mH + H2_traj * mH) / (mO + mH + mH)
        print("Calculating the gradient")
        self.com_water_velocity = np.gradient(self.com_water_traj, self.dtfs, axis=-1, edge_order=2)
        print("Calculating the autocorrelation function...")
        self.vacf_x, self.vacf_y, self.vacf_z = [], [], []
        for i in range(self.nwaters):
            self.vacf_x.append( self.auto_correlation_function_simple(self.com_water_velocity[i, 0, :]) )
            self.vacf_y.append( self.auto_correlation_function_simple(self.com_water_velocity[i, 1, :]) )
            self.vacf_z.append( self.auto_correlation_function_simple(self.com_water_velocity[i, 2, :]) )
        self.vacf_x = np.mean(np.array(self.vacf_x), axis=0)
        self.vacf_y = np.mean(np.array(self.vacf_y), axis=0)
        self.vacf_z = np.mean(np.array(self.vacf_z), axis=0)
        self.vacf_tot = self.vacf_x + self.vacf_y + self.vacf_z
        self.vacf_time_fs = np.linspace(0.0, self.dtfs*(self.vacf_x.size -1), self.vacf_x.size)
        print("Calculating the FFT of autocorrelation function")
        self.vacf_x_freq, self.vacf_x_sp = self.fft(self.vacf_x)
        self.vacf_y_freq, self.vacf_y_sp = self.fft(self.vacf_y)
        self.vacf_z_freq, self.vacf_z_sp = self.fft(self.vacf_z)
        self.vacf_tot_freq, self.vacf_tot_sp = self.fft(self.vacf_tot)

    def cacl_OHvelocity_water_traj(self):
        print("Calculating OH velocity traj of water")
        self.nwaters = self.natoms // 3
        mO, mH = 15.9994, 1.00794
        O_traj = self.traj[0:self.nwaters*3:3, :, :]
        H1_traj = self.traj[1:self.nwaters*3:3, :, :]
        H2_traj = self.traj[2:self.nwaters*3:3, :, :]
        self.OH_water_traj = O_traj - H1_traj
        print("Calculating the gradient")
        self.OH_water_velocity = np.gradient(self.OH_water_traj, self.dtfs, axis=-1, edge_order=2)
        print("Calculating the autocorrelation function...")
        self.vacf_OHx, self.vacf_OHy, self.vacf_OHz = [], [], []
        for i in range(self.nwaters):
            self.vacf_OHx.append( self.auto_correlation_function_simple(self.OH_water_velocity[i, 0, :]) )
            self.vacf_OHy.append( self.auto_correlation_function_simple(self.OH_water_velocity[i, 1, :]) )
            self.vacf_OHz.append( self.auto_correlation_function_simple(self.OH_water_velocity[i, 2, :]) )
        self.vacf_OHx = np.mean(np.array(self.vacf_OHx), axis=0)
        self.vacf_OHy = np.mean(np.array(self.vacf_OHy), axis=0)
        self.vacf_OHz = np.mean(np.array(self.vacf_OHz), axis=0)
        self.vacf_OHtot = self.vacf_OHx + self.vacf_OHy + self.vacf_OHz
        self.vacf_time_fs_OH = np.linspace(0.0, self.dtfs*(self.vacf_OHx.size -1), self.vacf_OHx.size)
        print("Calculating the FFT of autocorrelation function")
        self.vacf_OHx_freq, self.vacf_OHx_sp = self.fft(self.vacf_OHx)
        self.vacf_OHy_freq, self.vacf_OHy_sp = self.fft(self.vacf_OHy)
        self.vacf_OHz_freq, self.vacf_OHz_sp = self.fft(self.vacf_OHz)
        self.vacf_OHtot_freq, self.vacf_OHtot_sp = self.fft(self.vacf_OHtot)

    def cacl_dipole_water_traj(self, which_molecule=None, ntstart=0):
        print("Calculating dipole of water")
        cO, cH = -1.0, 0.5
        if which_molecule == None:
            self.nwaters = self.natoms // 3
            O_traj = self.traj[0:self.nwaters*3:3, :, ntstart:]
            H1_traj = self.traj[1:self.nwaters*3:3, :, ntstart:]
            H2_traj = self.traj[2:self.nwaters*3:3, :, ntstart:]
        else:
            O_traj = self.traj[3*which_molecule:3*which_molecule+3:3, :, ntstart:]
            H1_traj = self.traj[3*which_molecule+1:3*which_molecule+3:3, :, ntstart:]
            H2_traj = self.traj[3*which_molecule+2:3*which_molecule+3:3, :, ntstart:]
        self.dipole_water_traj = np.sum(O_traj * cO + H1_traj * cH + H2_traj * cH, axis=0)
        print("Calculating the autocorrelation function...")
        self.dacf_x = self.auto_correlation_function_simple(self.dipole_water_traj[0,:])
        self.dacf_y = self.auto_correlation_function_simple(self.dipole_water_traj[1,:])
        self.dacf_z = self.auto_correlation_function_simple(self.dipole_water_traj[2,:])
        self.dacf_tot = self.dacf_x + self.dacf_y + self.dacf_z
        self.dacf_time_fs = np.linspace(0.0, self.dtfs*(self.dacf_x.size -1), self.dacf_x.size)
        print("Calculating the FFT of autocorrelation function")
        self.dacf_x_freq, self.dacf_x_sp = self.fft3(self.dacf_x)
        self.dacf_y_freq, self.dacf_y_sp = self.fft3(self.dacf_y)
        self.dacf_z_freq, self.dacf_z_sp = self.fft3(self.dacf_z)
        self.dacf_tot_freq, self.dacf_tot_sp = self.fft3(self.dacf_tot)

    def cacl_dipole_water_traj_timespan(self, ntstart=0, ntend=10001, Nfft=10001):
        cO, cH = -1.0, 0.5
        self.nwaters = self.natoms // 3
        O_traj = self.traj[0:self.nwaters*3:3, :, ntstart:ntend]
        H1_traj = self.traj[1:self.nwaters*3:3, :, ntstart:ntend]
        H2_traj = self.traj[2:self.nwaters*3:3, :, ntstart:ntend]
        self.dipole_water_traj = np.sum(O_traj * cO + H1_traj * cH + H2_traj * cH, axis=0)
        self.dacf_x = self.auto_correlation_function_simple(self.dipole_water_traj[0,:])
        self.dacf_y = self.auto_correlation_function_simple(self.dipole_water_traj[1,:])
        self.dacf_z = self.auto_correlation_function_simple(self.dipole_water_traj[2,:])
        self.dacf_tot = self.dacf_x + self.dacf_y + self.dacf_z
        self.dacf_time_fs = np.linspace(0.0, self.dtfs*(self.dacf_x.size -1), self.dacf_x.size)
        self.dacf_x_freq, self.dacf_x_sp = self.fft3(self.dacf_x, N=Nfft)
        self.dacf_y_freq, self.dacf_y_sp = self.fft3(self.dacf_y, N=Nfft)
        #self.dacf_z_freq, self.dacf_z_sp = self.fft3(self.dacf_z, N=Nfft)
        #self.dacf_tot_freq, self.dacf_tot_sp = self.fft3(self.dacf_tot, N=Nfft)
        return self.dacf_x_freq, self.dacf_x_sp, self.dacf_y_sp, self.dacf_y_freq

    def calc_dipole_timeresolved(self, ntstart=0, ntend=10001):
        cO, cH = -1.0, 0.5
        self.nwaters = self.natoms // 3
        O_traj = self.traj[0:self.nwaters*3:3, :, ntstart:ntend]
        H1_traj = self.traj[1:self.nwaters*3:3, :, ntstart:ntend]
        H2_traj = self.traj[2:self.nwaters*3:3, :, ntstart:ntend]
        dipole_water_traj = np.sum(O_traj * cO + H1_traj * cH + H2_traj * cH, axis=0)
        dx, dy, dz = dipole_water_traj[0,:], dipole_water_traj[1,:], dipole_water_traj[2,:]
        dtime_fs = np.linspace(0.0, self.dtfs*(dx.size -1),  dx.size)
        return dtime_fs, dx, dy, dz

    def cacl_photon_traj_timespan(self, ntstart=0, ntend=-1):
        ph1 = self.traj[-2, 0, ntstart:ntend]
        ph2 = self.traj[-1, 1, ntstart:ntend]
        ph_time_fs = np.linspace(0.0, self.dtfs*(ph1.size -1), ph1.size)
        # I also want to calculate the velocity of these two photons
        v_ph1 = np.gradient(ph1, self.dtfs, axis=-1, edge_order=2)
        v_ph2 = np.gradient(ph2, self.dtfs, axis=-1, edge_order=2)
        return ph_time_fs, ph1, ph2, v_ph1, v_ph2
    
    def cacl_local_IR_traj_timespan(self, whichmolecule=0, ntstart=0, ntend=10001, Nfft=10001):
        cO, cH = -1.0, 0.5
        O_traj = self.traj[whichmolecule*3, :, ntstart:ntend]
        H1_traj = self.traj[whichmolecule*3+1, :, ntstart:ntend]
        H2_traj = self.traj[whichmolecule*3+2, :, ntstart:ntend]
        self.dipole_water_traj = O_traj * cO + H1_traj * cH + H2_traj * cH
        self.dacf_x = self.auto_correlation_function_simple(self.dipole_water_traj[0,:])
        self.dacf_y = self.auto_correlation_function_simple(self.dipole_water_traj[1,:])
        self.dacf_z = self.auto_correlation_function_simple(self.dipole_water_traj[2,:])
        self.dacf_tot = self.dacf_x + self.dacf_y + self.dacf_z
        self.dacf_time_fs = np.linspace(0.0, self.dtfs*(self.dacf_x.size -1), self.dacf_x.size)
        self.dacf_x_freq, self.dacf_x_sp = self.fft3(self.dacf_x, N=Nfft)
        self.dacf_y_freq, self.dacf_y_sp = self.fft3(self.dacf_y, N=Nfft)
        self.dacf_z_freq, self.dacf_z_sp = self.fft3(self.dacf_z, N=Nfft)
        #self.dacf_tot_freq, self.dacf_tot_sp = self.fft3(self.dacf_tot, N=Nfft)
        return self.dacf_x_freq, self.dacf_x_sp, self.dacf_y_sp, self.dacf_z_sp

    def cacl_CO_energy_traj_timespan(self, ntstart=0, ntend=10001, nskipping=1, nmol=216, summed=False):
        C_traj = self.traj[0:nmol*3:3, :, ntstart:ntend]
        O1_traj = self.traj[1:nmol*3:3, :, ntstart:ntend]
        O2_traj = self.traj[2:nmol*3:3, :, ntstart:ntend]
        CO_1 = np.sqrt(np.sum(np.abs(C_traj - O1_traj)**2, axis=1)) * 1.8897259885789 # to a.u.
        CO_2 = np.sqrt(np.sum(np.abs(C_traj - O2_traj)**2, axis=1)) * 1.8897259885789 # to a.u.
        def VCO(x):
            x -= 2.196
            return 0.45106005 * x**2 - 0.67297933 * x**3 + 0.5857145 * x**4
        VCO_smoothed = VCO(CO_1) + VCO(CO_2)
        print("VCO shape", np.shape(VCO_smoothed))
        for i in range(VCO_smoothed.shape[0]):
            VCO_smoothed[i,:] = smooth(VCO_smoothed[i,:])
        time_fs = np.linspace(0.0, self.dtfs*(VCO_smoothed.shape[-1] -1), VCO_smoothed.shape[-1])[::nskipping]
        VCO_smoothed = VCO_smoothed[:, ::nskipping]
        print("VCO shape", np.shape(VCO_smoothed))
        if summed:
            VCO_smoothed = np.reshape(np.sum(VCO_smoothed, axis=0), (1, -1))
        print("VCO shape", np.shape(VCO_smoothed))
        return time_fs, VCO_smoothed


    def cacl_orientation_water_traj(self):
        print("Calculating orientation correlation function of water")
        self.nwaters = self.natoms // 3
        #self.cacl_center_of_mass_water_traj()
        mO, mH = 15.9994, 1.00794
        O_traj = self.traj[0:self.nwaters*3:3, :, :]
        H1_traj = self.traj[1:self.nwaters*3:3, :, :]
        H2_traj = self.traj[2:self.nwaters*3:3, :, :]
        # Calculate the unit vector along the dipole moment direction
        self.dipoles = -O_traj * 2.0 + H1_traj + H2_traj
        self.dipoles_length = np.sqrt(np.sum(self.dipoles**2, axis=1))
        self.dipoles[:,0,:] /= self.dipoles_length
        self.dipoles[:,1,:] /= self.dipoles_length
        self.dipoles[:,2,:] /= self.dipoles_length
        print("Calculating the autocorrelation function...")
        self.oacfz_x, self.oacfz_y, self.oacfz_z = [], [], []
        for i in range(self.nwaters):
            self.oacfz_x.append( self.auto_correlation_function_simple(self.dipoles[i, 0, :]) )
            self.oacfz_y.append( self.auto_correlation_function_simple(self.dipoles[i, 1, :]) )
            self.oacfz_z.append( self.auto_correlation_function_simple(self.dipoles[i, 2, :]) )
        self.oacfz_x = np.mean(np.array(self.oacfz_x), axis=0)
        self.oacfz_y = np.mean(np.array(self.oacfz_y), axis=0)
        self.oacfz_z = np.mean(np.array(self.oacfz_z), axis=0)
        self.oacfz_tot = self.oacfz_x + self.oacfz_y + self.oacfz_z
        self.oacfz_time_fs = np.linspace(0.0, self.dtfs*(self.oacfz_x.size -1), self.oacfz_x.size)
        print("Calculating the FFT of autocorrelation function")
        self.oacfz_x_freq, self.oacfz_x_sp = self.fft3(self.oacfz_x)
        self.oacfz_y_freq, self.oacfz_y_sp = self.fft3(self.oacfz_y)
        self.oacfz_z_freq, self.oacfz_z_sp = self.fft3(self.oacfz_z)
        self.oacfz_tot_freq, self.oacfz_tot_sp = self.fft3(self.oacfz_tot)

    def calc_bond_length_dist(self, bin_start=0.1, bin_end=4, bin_num=1000, ntstart=0, ntend=-1):
        self.nwaters = self.natoms // 3
        O_traj = self.traj[0:self.nwaters*3:3, :, ntstart:ntend]
        H1_traj = self.traj[1:self.nwaters*3:3, :, ntstart:ntend]
        H2_traj = self.traj[2:self.nwaters*3:3, :, ntstart:ntend]
        bond_bins = np.linspace(bin_start, bin_end, bin_num)
        print("Calculating Bond length distribution of water")
        data1 = np.reshape(np.sqrt(np.sum(np.abs(O_traj - H1_traj)**2, axis=1)), -1)
        data2 = np.reshape(np.sqrt(np.sum(np.abs(O_traj - H2_traj)**2, axis=1)), -1)
        data = np.concatenate((data1, data2))
        print("preparing data")
        self.bond_statistics, self.bond_bins = np.histogram(data, bins=bond_bins)
        return self.bond_bins, self.bond_statistics

    def calc_bond_length_fluc(self, ntstart=0, ntend=-1):
        self.nwaters = self.natoms // 3
        O_traj = self.traj[0:self.nwaters*3:3, :, ntstart:ntend]
        H1_traj = self.traj[1:self.nwaters*3:3, :, ntstart:ntend]
        H2_traj = self.traj[2:self.nwaters*3:3, :, ntstart:ntend]
        data1 = np.sqrt(np.sum(np.abs(O_traj - H1_traj)**2, axis=1))
        data2 = np.sqrt(np.sum(np.abs(O_traj - H2_traj)**2, axis=1))
        data = np.concatenate((data1, data2), axis=0)
        # I have calculated the bond length of C=O in the time series
        data_std = np.std(data, axis=0)
        time_fs = np.linspace(0.0, self.dtfs*(data_std.size -1), data_std.size)
        return time_fs, smooth(data_std)
    
    def calc_pair_dist_OO(self, bin_start=0.1, bin_end=9.0, bin_num=1000, cell_length=18.6445):
        self.nwaters = self.natoms // 3
        O_traj = self.traj[0:self.nwaters*3:3, :, 1:-1:10]
        H1_traj = self.traj[1:self.nwaters*3:3, :, 1:-1:10]
        H2_traj = self.traj[2:self.nwaters*3:3, :, 1:-1:10]
        r_bins = np.linspace(bin_start, bin_end, bin_num)
        print("Calculating pari distribution function of water O-O")
        self.OOstatstics = np.zeros(np.size(r_bins)-1)
        for i in range(self.nwaters):
            for j in range(self.nwaters):
                OO_traj = np.abs(O_traj[i, :, :] - O_traj[j, :, :]) % cell_length
                OO_traj[OO_traj >= cell_length / 2.0] -= cell_length
                distance_oo = np.sqrt(np.sum(OO_traj**2, axis=0))
                oo, self.r_bins = np.histogram(distance_oo, bins=r_bins)
                oo = oo.astype(np.float)
                self.OOstatstics += oo
        # post process
        self.OOstatstics /= float(self.nwaters)
        self.OOstatstics /= 4.0 * np.pi * r_bins[1:]**2 * (r_bins[2] - r_bins[1])
        self.OOstatstics /= float(self.nwaters) / cell_length**3
        self.OOstatstics /= np.shape(O_traj)[-1]

    def auto_correlation_function_fft(self, x):
        corr = signal.fftconvolve(x, x[::-1], mode='same')
        corr = corr[corr.size // 2: ]
        return corr / corr[0] * np.mean(x * x)

    def auto_correlation_function_simple(self, x):
        n = x.size
        if n % 2 == 0:
            x_shifted = np.zeros(n*2)
        else:
            x_shifted = np.zeros(n*2-1)
        x_shifted[n//2 : n//2+n] = x
        # Convolute the shifted array with the flipped array, which is equivalent to performing a correlation
        autocorr_full = (signal.fftconvolve(x_shifted, x[::-1], mode='same')[-n:]/ np.arange(n, 0, -1))
        # Truncate the autocorrelation array
        autocorr = autocorr_full[0:n//2]
        return autocorr

    def fft(self, x):
        #sp = np.fft.fft(x)
        #freq_au = 2.0 * np.pi * np.fft.fftfreq(np.size(x), self.dtau)
        # Because dt has the unit of fs, I need to transform fs^{-1} to cm^{-1}
        #freq_cminverse = freq_au * 219474.63
        #return freq_cminverse[0:sp.size//2], sp[0:sp.size//2]
        lineshape = fftpack.dct(x, type=1)
        freq_au = np.linspace(0, 0.5/self.dtfs * 1e15, len(x))
        # Because dt has the unit of fs, I need to transform fs^{-1} to cm^{-1}
        freq_cminverse = freq_au / (100.0 * 299792458.0)
        return freq_cminverse, lineshape

    def fft3(self, x, N=None):
        # Adding zeros to the end of x
        if N is not None:
            n = N
        else:
            n = np.size(x)
        lineshape = fftpack.dct(x, type=1, n=n)
        freq_au = np.linspace(0, 0.5/self.dtfs * 1e15, n)
        # Because dt has the unit of fs, I need to transform fs^{-1} to cm^{-1}
        freq_cminverse = freq_au / (100.0 * 299792458.0)
        # Calculate spectra
        #field_description =  freq_au**2
        field_description =  freq_au**2
        spectra = lineshape * field_description
        return freq_cminverse, spectra
        #return freq_cminverse[0:spectra.size//2], spectra[0:spectra.size//2]

    def output_velocity_autocorrelation(self):
        local_filename = "%s.vac.txt" %self.xyz_filename
        if os.path.isfile(local_filename + "dd"):
            print("Have calculated COM diffusion for %s, skipping..." %self.xyz_filename)
        else:
            self.cacl_center_of_mass_water_traj()
            # output data
            data = np.zeros((np.size(self.vacf_time_fs), 10))
            data[:, 0] = self.vacf_time_fs
            data[:, 1] = self.vacf_x
            data[:, 2] = self.vacf_y
            data[:, 3] = self.vacf_z
            data[:, 4] = self.vacf_tot
            data[:, 5] = self.vacf_x_freq
            data[:, 6] = smooth(self.vacf_x_sp)
            data[:, 7] = smooth(self.vacf_y_sp)
            data[:, 8] = smooth(self.vacf_z_sp)
            data[:, 9] = smooth(self.vacf_tot_sp)
            comments = "# vacf_time_fs, vacf_x, vacf_y, vacf_z, vacf_tot, freq, sp_x, sp_y, sp_z, sp_tot"
            np.savetxt(local_filename, data, comments=comments)

    def output_OHvelocity_autocorrelation(self):
        local_filename = "%s.OHvac.txt" %self.xyz_filename
        if os.path.isfile(local_filename):
            print("Have calculated OH bond diffusion for %s, skipping..." %self.xyz_filename)
        else:
            self.cacl_OHvelocity_water_traj()
            # output data
            data = np.zeros((np.size(self.vacf_time_fs_OH), 10))
            data[:, 0] = self.vacf_time_fs_OH
            data[:, 1] = self.vacf_OHx
            data[:, 2] = self.vacf_OHy
            data[:, 3] = self.vacf_OHz
            data[:, 4] = self.vacf_OHtot
            data[:, 5] = self.vacf_OHx_freq
            data[:, 6] = smooth(self.vacf_OHx_sp)
            data[:, 7] = smooth(self.vacf_OHy_sp)
            data[:, 8] = smooth(self.vacf_OHz_sp)
            data[:, 9] = smooth(self.vacf_OHtot_sp)
            comments = "# vacf_time_fs, vacf_x, vacf_y, vacf_z, vacf_tot, freq, sp_x, sp_y, sp_z, sp_tot"
            np.savetxt(local_filename, data, comments=comments)

    def output_dipole_autocorrelation(self, which_molecule=None, ntstart=0):
        if which_molecule is None:
            local_filename = "%s.dac.txt" %self.xyz_filename
        else:
            local_filename = "%s.dac_%d.txt" %(self.xyz_filename, which_molecule)
        if os.path.isfile(local_filename):
            print("Have calculated dipole autocorrelation for %s, skipping..." %self.xyz_filename)
        else:
            self.cacl_dipole_water_traj(which_molecule=which_molecule, ntstart=ntstart)
            # output data
            data = np.zeros((np.size(self.dacf_x), 10))
            data[:, 0] = self.dacf_time_fs
            data[:, 1] = self.dacf_x
            data[:, 2] = self.dacf_y
            data[:, 3] = self.dacf_z
            data[:, 4] = self.dacf_tot
            data[:, 5] = self.dacf_x_freq
            data[:, 6] = smooth(self.dacf_x_sp)
            data[:, 7] = smooth(self.dacf_y_sp)
            data[:, 8] = smooth(self.dacf_z_sp)
            data[:, 9] = smooth(self.dacf_tot_sp)
            comments = "dacf_time_fs, dacf_x, dacf_y, dacf_z, dacf_tot, freq, sp_x, sp_y, sp_z, sp_tot"
            np.savetxt(local_filename, data, comments=comments, fmt='%.4E')

    def output_dipole_time_resolved(self, ntstart=0, ntend=10001):
        local_filename = "%s.dipole_td.txt" %self.xyz_filename
        if os.path.isfile(local_filename + "dd"):
            print("Have calculated dipole time resolved for %s, skipping..." %self.xyz_filename)
        else:
            t, dx, dy, dz = self.calc_dipole_timeresolved(ntstart=ntstart, ntend=ntend)
            # output data
            data = np.zeros((np.size(dx), 4))
            data[:, 0] = t
            data[:, 1] = dx
            data[:, 2] = dy
            data[:, 3] = dz
            comments = "time fs, dx, dy, dz"
            np.savetxt(local_filename, data, comments=comments, fmt="%.4E")

    def output_IR_dynamics(self, nspan=10, nresolution=10, nfft=2501):
        local_filename = "%s.dacf_td.txt" %self.xyz_filename
        if os.path.isfile(local_filename):
            print("Have calculated time-dependent dipole autocorrelation for %s, skipping..." %self.xyz_filename)
        else:
            print("performing time-dependent IR calculation...")
            # Firstly, try one
            N = 50001 // nspan
            dN = 50001 // nresolution
            freq, spx, spy, spz = self.cacl_dipole_water_traj_timespan(ntstart=0, ntend=N, Nfft=nfft)
            data = np.zeros((np.size(freq), 1+nresolution))
            data[:, 0] = freq
            for i in range(nresolution):
                freq, spx, spy, spz = self.cacl_dipole_water_traj_timespan(ntstart=dN*i, ntend=dN*i+N, Nfft=nfft)
                data[:, 1+i] = smooth(spx + spy)
            np.savetxt(local_filename, data, fmt='%.4E')

    def output_local_IR_dynamics(self, nspan=10, nresolution=10, nfft=10001, nmol_start=0, nmol_end=216):
        local_filename = "%s.localIR_td.txt" %self.xyz_filename
        if nmol_end - nmol_start is not 216:
            local_filename = "%s.localIR_%d_%d_td.txt" %(self.xyz_filename, nmol_start, nmol_end)
        if os.path.isfile(local_filename):
            print("Have calculated time-dependent local IR for %s, skipping..." %self.xyz_filename)
        else:
            # Firstly, try one
            N = 50001 // nspan
            dN = 50001 // nresolution
            freq, spx, spy, spz = self.cacl_local_IR_traj_timespan(whichmolecule=0, ntstart=0, ntend=N, Nfft=nfft)
            data = np.zeros((np.size(freq), 1+nresolution))
            # After obtaining the result for the first molecule, we repeat the procedure to calculate the remaining molecules
            self.nwaters = 216
            if nmol_end > self.nwaters:
                nmol_end = self.nwaters
            for n in range(nmol_start, nmol_end):
                if n%10 == 0:
                    print("performing time-dependent local IR calculation for No.%d water when %d" %(n, nmol_end-nmol_start))
                for i in range(nresolution):
                    freq, spx, spy, spz = self.cacl_local_IR_traj_timespan(whichmolecule=n, ntstart=dN*i, ntend=dN*i+N, Nfft=nfft)
                    data[:, 1+i] += smooth(spx + spy + spz)
                data[:,0] += freq
            data /= (nmol_end - nmol_start)
            np.savetxt(local_filename, data, fmt='%.4E')

    def output_ph_energy_dynamics(self, ntstart=0, ntend=50001):
        local_filename = "%s.ph_td.txt" %self.xyz_filename
        if os.path.isfile(local_filename):
            print("Have calculated dynamics of photons for %s, skipping..." %self.xyz_filename)
        else:
            t, ph1, ph2, v_ph1, v_ph2 = self.cacl_photon_traj_timespan(ntstart=ntstart, ntend=ntend)
            t, bond_std = self.calc_bond_length_fluc(ntstart=ntstart, ntend=ntend)
            data = np.zeros((np.size(t), 6))
            data[:,0] = t
            data[:,1] = ph1
            data[:,2] = ph2
            data[:,3] = v_ph1
            data[:,4] = v_ph2
            data[:,5] = bond_std
            np.savetxt(local_filename, data, fmt='%.4E')
    
    def output_bond_length_dynamics(self, bin_start=0.5, bin_end=1.5, bin_num=1000, nspan=10, nresolution=10):
        local_filename = "%s.bond_length_dist_td.txt" %self.xyz_filename
        if os.path.isfile(local_filename + "dd"):
            print("Have calculated bond length distribution dynamics for %s, skipping..." %self.xyz_filename)
        else:
            N = 50001 // nspan
            dN = 50001 // nresolution
            bins, count = self.calc_bond_length_dist(bin_start=bin_start, bin_end=bin_end, bin_num=bin_num, ntstart=0, ntend=N)
            data = np.zeros((np.size(bins)-1, 1+nresolution))
            # output data
            for i in range(nresolution):
                bins, count = self.calc_bond_length_dist(bin_start=bin_start, bin_end=bin_end, bin_num=bin_num, ntstart=dN*i, ntend=dN*i+N)
                data[:, 1+i] = count.astype(np.float) #/ np.sum(count)
            data[:, 0] = bins[:-1]
            comments = "# bins, counts"
            np.savetxt(local_filename, data, comments=comments, fmt='%.4E')

    def output_VCO_dynamics(self, ntstart=0, ntend=10001, nskipping=1, summed=False, nmol=216):
        local_filename = "%s.VCO_td.txt" %self.xyz_filename
        if summed:
            local_filename = "%s.VCO_summed_td.txt" %self.xyz_filename
        if os.path.isfile(local_filename + "dd"):
            print("Have calculated VCO dynamics for %s, skipping..." %self.xyz_filename)
        else:
            time_fs, VCO_traj = self.cacl_CO_energy_traj_timespan(ntstart=ntstart, ntend=ntend, nskipping=nskipping, summed=summed, nmol=nmol)
            data = np.zeros((VCO_traj.shape[1], VCO_traj.shape[0]+1))
            data[:,0] = time_fs
            data[:,1:] = np.transpose(VCO_traj)
            comments = "# time, VCO for each molecule"
            np.savetxt(local_filename, data, comments=comments, fmt='%.4E')
    
    
    def output_local_IR_statistics(self, ntstart=0, ntend=2501, nfft=2501):
        local_filename = "%s.localIR_statistics_%d_%d.txt" %(self.xyz_filename, ntstart, ntend)
        if os.path.isfile(local_filename):
            print("Have calculated local IR statistics for %s, skipping..." %self.xyz_filename)
        self.nwaters = 216
        freq, spx, spy, spz = self.cacl_local_IR_traj_timespan(whichmolecule=0, ntstart=ntstart, ntend=ntend, Nfft=nfft)
        data = np.zeros((np.size(freq), 1 + self.nwaters))
        for n in range(self.nwaters):
            freq, spx, spy, spz = self.cacl_local_IR_traj_timespan(whichmolecule=n, ntstart=ntstart, ntend=ntend, Nfft=nfft)
            data[:,1+n] = smooth(spx + spy + spz)
        data[:,0] = freq
        np.savetxt(local_filename, data, fmt='%.4E')

    def output_orientation_autocorrelation(self):
        local_filename = "%s.oac1.txt" %self.xyz_filename
        if os.path.isfile(local_filename + "dd"):
            print("Have calculated 1st orientation autocorrelation for %s, skipping..." %self.xyz_filename)
        else:
            self.cacl_orientation_water_traj()
            # output data
            data = np.zeros((np.size(self.oacfz_time_fs), 10))
            data[:, 0] = self.oacfz_time_fs
            data[:, 1] = self.oacfz_x
            data[:, 2] = self.oacfz_y
            data[:, 3] = self.oacfz_z
            data[:, 4] = self.oacfz_tot
            data[:, 5] = self.oacfz_x_freq
            data[:, 6] = smooth(self.oacfz_x_sp)
            data[:, 7] = smooth(self.oacfz_y_sp)
            data[:, 8] = smooth(self.oacfz_z_sp)
            data[:, 9] = smooth(self.oacfz_tot_sp)
            comments = "# oacfz_time_fs, oacfz_x, oacfz_y, oacfz_z, oacfz_tot, freq, sp_x, sp_y, sp_z, sp_tot"
            np.savetxt(local_filename, data, comments=comments)

    def output_bond_length(self, bin_start=0.5, bin_end=1.5, bin_num=1000):
        local_filename = "%s.bond_length_dist.txt" %self.xyz_filename
        if os.path.isfile(local_filename + "dd"):
            print("Have calculated bond length distribution for %s, skipping..." %self.xyz_filename)
        else:
            self.calc_bond_length_dist(bin_start=bin_start, bin_end=bin_end, bin_num=bin_num)
            # output data
            data = np.zeros((np.size(self.bond_bins)-1, 2))
            data[:, 0] = self.bond_bins[:-1]
            data[:, 1] = self.bond_statistics.astype(np.float)
            comments = "# bins, counts"
            np.savetxt(local_filename, data, comments=comments)

    def output_pair_dist(self, bin_start=0.1, bin_end=9, bin_num=1000, cell_length=18.6445):
        local_filename = "%s.pair_dist.txt" %self.xyz_filename
        if os.path.isfile(local_filename + "dd"):
            print("Have calculated pair distribution function for %s, skipping..." %self.xyz_filename)
        else:
            self.calc_pair_dist_OO(bin_start=bin_start, bin_end=bin_end, bin_num=bin_num, cell_length=cell_length)
            # output data
            data = np.zeros((np.size(self.r_bins)-1, 2))
            data[:, 0] = self.r_bins[1:]
            data[:, 1] = self.OOstatstics.astype(np.float)
            comments = "# bins, OO"
            print("output data as %s" %local_filename)
            np.savetxt(local_filename, data, comments=comments)


if __name__ == "__main__":
    command = sys.argv[1]
    if command == "-dacf":
        default_file_size=313681365
    elif command == "-mu":
        default_file_size=313681365
    elif command == "-single":
        default_file_size=4030403
    else:
        default_file_size=1568281365
    paths = sys.argv[2:]
    for i, path in enumerate(paths):
        filenames = glob.glob("%s/simu*.xc.xyz" %path)
        for j, filename in enumerate(filenames):
            if os.path.getsize(filename) >= default_file_size:
                if command == "-dacf"  or command == "-single":
                    a = MD_Analysis(xyz_filename=filename, nframe_max=10001)
                else:
                    a = MD_Analysis(xyz_filename=filename, nframe_max=50001)
                if command == "-mu":
                    a.output_dipole_time_resolved()
                if command == "-dacf":
                    a.output_dipole_autocorrelation(which_molecule=None)
                if command == "-single":
                    a.output_dipole_autocorrelation(which_molecule=None)
                if command == "-ph" or command == "-phIR":
                    a.output_ph_energy_dynamics()
                if command == "-vco":
                    a.output_VCO_dynamics(ntstart=0, ntend=1001, nskipping=25)
                    #a.output_VCO_dynamics(ntstart=0, ntend=50001, nskipping=200)
                    #a.output_VCO_dynamics(ntstart=0, ntend=50001, nskipping=10, summed=True)
                if command == "-vco1":
                    a.output_VCO_dynamics(nmol=1, ntend=50001)
                #a.output_VCO_dynamics(ntstart=0, ntend=2501, nskipping=1)
                #a.output_VCO_dynamics(ntstart=0, ntend=30000, nskipping=1, summed=True)
                #a.output_VCO_dynamics(ntstart=0, ntend=50001, nskipping=100)
                #if j == 0:
                #    for n in [2500*i for i in range(10)]:
                #    for n in [100*i for i in range(50)]:
                #        a.output_local_IR_statistics(ntstart=n, ntend=2501+n, nfft=2501)
                #    n = 500
                #    a.output_local_IR_statistics(ntstart=n, ntend=2501+n, nfft=2501)
                if command == "-IR" or command == "-phIR":
                    a.output_IR_dynamics(nspan=20, nresolution=100, nfft=2501)
                    a.output_local_IR_dynamics(nspan=20, nresolution=100, nfft=2501)
                #a.output_local_IR_dynamics(nspan=20, nresolution=100, nfft=2501, nmol_start=0, nmol_end=107)
                #a.output_local_IR_dynamics(nspan=20, nresolution=100, nfft=2501, nmol_start=107, nmol_end=216)
