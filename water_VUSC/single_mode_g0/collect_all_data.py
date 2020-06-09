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

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."

    if window_len<3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[window_len//2-1:-window_len//2]

class MD_Analysis:
    def __init__(self, xyz_filename, dtfs = 2, nframe_max=10001):
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
        mO, mH = 15.9994, 1.00794
        O_traj = self.traj[0:self.nwaters*3:3, :, :]
        H1_traj = self.traj[1:self.nwaters*3:3, :, :]
        H2_traj = self.traj[2:self.nwaters*3:3, :, :]
        self.com_water_traj = (O_traj * mO + H1_traj * mH + H2_traj * mH) / (mO + 2.0 * mH)
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

    def cacl_Ovelocity_water_traj(self):
        print("Calculating Oxygen velocity traj of water")
        self.nwaters = self.natoms // 3
        mO, mH = 15.9994, 1.00794
        O_traj = self.traj[0:self.nwaters*3:3, :, :]
        print("Calculating the gradient")
        self.O_water_velocity = np.gradient(O_traj, self.dtfs, axis=-1, edge_order=2)
        print("Calculating the autocorrelation function...")
        self.vacf_Ox, self.vacf_Oy, self.vacf_Oz = [], [], []
        for i in range(self.nwaters):
            self.vacf_Ox.append( self.auto_correlation_function_simple(self.O_water_velocity[i, 0, :]) )
            self.vacf_Oy.append( self.auto_correlation_function_simple(self.O_water_velocity[i, 1, :]) )
            self.vacf_Oz.append( self.auto_correlation_function_simple(self.O_water_velocity[i, 2, :]) )
        self.vacf_Ox = np.mean(np.array(self.vacf_Ox), axis=0)
        self.vacf_Oy = np.mean(np.array(self.vacf_Oy), axis=0)
        self.vacf_Oz = np.mean(np.array(self.vacf_Oz), axis=0)
        self.vacf_Otot = self.vacf_Ox + self.vacf_Oy + self.vacf_Oz
        self.vacf_time_fs_O = np.linspace(0.0, self.dtfs*(self.vacf_Ox.size -1), self.vacf_Ox.size)
        print("Calculating the FFT of autocorrelation function")
        self.vacf_Ox_freq, self.vacf_Ox_sp = self.fft(self.vacf_Ox)
        self.vacf_Oy_freq, self.vacf_Oy_sp = self.fft(self.vacf_Oy)
        self.vacf_Oz_freq, self.vacf_Oz_sp = self.fft(self.vacf_Oz)
        self.vacf_Otot_freq, self.vacf_Otot_sp = self.fft(self.vacf_Otot)

    def cacl_Hvelocity_water_traj(self):
        print("Calculating Hydrogen velocity traj of water")
        self.nwaters = self.natoms // 3
        mO, mH = 15.9994, 1.00794
        H1_traj = self.traj[1:self.nwaters*3:3, :, :]
        print("Calculating the gradient")
        self.H_water_velocity = np.gradient(H1_traj, self.dtfs, axis=-1, edge_order=2)
        print("Calculating the autocorrelation function...")
        self.vacf_Hx, self.vacf_Hy, self.vacf_Hz = [], [], []
        for i in range(self.nwaters):
            self.vacf_Hx.append( self.auto_correlation_function_simple(self.H_water_velocity[i, 0, :]) )
            self.vacf_Hy.append( self.auto_correlation_function_simple(self.H_water_velocity[i, 1, :]) )
            self.vacf_Hz.append( self.auto_correlation_function_simple(self.H_water_velocity[i, 2, :]) )
        self.vacf_Hx = np.mean(np.array(self.vacf_Hx), axis=0)
        self.vacf_Hy = np.mean(np.array(self.vacf_Hy), axis=0)
        self.vacf_Hz = np.mean(np.array(self.vacf_Hz), axis=0)
        self.vacf_Htot = self.vacf_Hx + self.vacf_Hy + self.vacf_Hz
        self.vacf_time_fs_H = np.linspace(0.0, self.dtfs*(self.vacf_Hx.size -1), self.vacf_Hx.size)
        print("Calculating the FFT of autocorrelation function")
        self.vacf_Hx_freq, self.vacf_Hx_sp = self.fft(self.vacf_Hx)
        self.vacf_Hy_freq, self.vacf_Hy_sp = self.fft(self.vacf_Hy)
        self.vacf_Hz_freq, self.vacf_Hz_sp = self.fft(self.vacf_Hz)
        self.vacf_Htot_freq, self.vacf_Htot_sp = self.fft(self.vacf_Htot)

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

    def cacl_OHkinetic_water_traj(self):
        print("Calculating OH velocity traj of water")
        self.nwaters = self.natoms // 3
        mO, mH = 15.9994, 1.00794
        cO, cH = -1.0, 0.5
        O_traj = self.traj[0:self.nwaters*3:3, :, :]
        H1_traj = self.traj[1:self.nwaters*3:3, :, :]
        H2_traj = self.traj[2:self.nwaters*3:3, :, :]
        print("Calculating the gradient")
        # Do coordinate transform
        OH1 = H1_traj - O_traj
        OH2 = H2_traj - O_traj
        OH1_length = np.sqrt(np.sum(OH1**2, axis=1))
        OH2_length = np.sqrt(np.sum(OH2**2, axis=1))
        theta = np.arccos(np.sum(OH1 * OH2, axis=1)) / OH1_length / OH2_length / 2.0
        
        H1 = np.zeros(np.shape(O_traj))
        H2 = np.zeros(np.shape(O_traj))
        H1[:,1,:] = OH1_length * np.sin(theta)
        H1[:,2,:] = OH1_length * np.cos(theta)
        H2[:,1,:] = OH2_length * np.sin(-theta)
        H2[:,2,:] = OH2_length * np.cos(-theta)

        self.OH_water_traj = np.zeros(np.shape(O_traj))
        self.OH_water_traj[:, 1, :] = H1[:, 1, :] + H2[:, 1, :]
        self.OH_water_traj[:, 2, :] = H1[:, 2, :] + H2[:, 2, :]

        self.OH_water_velocity = np.gradient(self.OH_water_traj, self.dtfs, axis=-1, edge_order=2)**2
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

    def cacl_dipole_water_traj(self):
        print("Calculating dipole of water")
        self.nwaters = self.natoms // 3
        cO, cH = -1.0, 0.5
        O_traj = self.traj[0:self.nwaters*3:3, :, :]
        H1_traj = self.traj[1:self.nwaters*3:3, :, :]
        H2_traj = self.traj[2:self.nwaters*3:3, :, :]
        self.dipole_water_traj = np.sum(O_traj * cO + H1_traj * cH + H2_traj * cH, axis=0)
        #self.dipole_water_traj = np.gradient(self.dipole_water_traj, self.dtfs, axis=-1, edge_order=2)
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

    def cacl_orientation_water_traj2(self):
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
        xx, xy, xz, yy, yz, zz = [], [], [], [], [], []
        for i in range(self.nwaters):
            xx.append( self.auto_correlation_function_simple(self.dipoles[i, 0, :]**2) )
            xy.append( self.auto_correlation_function_simple(self.dipoles[i, 0, :] * self.dipoles[i, 1, :]) )
            xz.append( self.auto_correlation_function_simple(self.dipoles[i, 0, :] * self.dipoles[i, 2, :]) )
            yy.append( self.auto_correlation_function_simple(self.dipoles[i, 1, :]**2) )
            yz.append( self.auto_correlation_function_simple(self.dipoles[i, 1, :] * self.dipoles[i, 2, :]) )
            zz.append( self.auto_correlation_function_simple(self.dipoles[i, 1, :]**2) )
        xx = np.mean(np.array(xx), axis=0)
        xy = np.mean(np.array(xy), axis=0)
        xz = np.mean(np.array(xz), axis=0)
        yy = np.mean(np.array(yy), axis=0)
        yz = np.mean(np.array(yz), axis=0)
        zz = np.mean(np.array(zz), axis=0)
        self.oacfz2_tot = 1.5*(xx + xy*2.0 + xz*2.0 + yy + yz*2.0 + zz) - 0.5
        self.oacfz2_time_fs = np.linspace(0.0, self.dtfs*(self.oacfz2_tot.size -1), self.oacfz2_tot.size)
        print("Calculating the FFT of autocorrelation function")
        self.oacfz2_tot_freq, self.oacfz2_tot_sp = self.fft3(self.oacfz2_tot)

    def calc_bond_length_dist(self, bin_start=0.1, bin_end=4, bin_num=1000):
        self.nwaters = self.natoms // 3
        O_traj = self.traj[0:self.nwaters*3:3, :, :]
        H1_traj = self.traj[1:self.nwaters*3:3, :, :]
        H2_traj = self.traj[2:self.nwaters*3:3, :, :]
        bond_bins = np.linspace(bin_start, bin_end, bin_num)
        print("Calculating Bond length distribution of water")
        OH = O_traj - H1_traj
        data = np.sqrt(np.sum(np.abs(O_traj - H1_traj)**2, axis=1))
        print("preparing data")
        self.bond_statistics, self.bond_bins = np.histogram(data, bins=bond_bins)

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

    def fft2(self, x ):
        # Adding zeros to the end of x
        lineshape = fftpack.dct(x, type=1)
        freq_au = np.linspace(0, 0.5/self.dtfs * 1e15, len(x))
        # Because dt has the unit of fs, I need to transform fs^{-1} to cm^{-1}
        freq_cminverse = freq_au / (100.0 * 299792458.0)
        # Calculate spectra
        #field_description =  freq_au**2
        KT = 0.0009500372557109825
        field_description =  freq_au*(1.0 - np.exp(-freq_au / KT))
        spectra = lineshape * field_description
        return freq_cminverse, spectra
        #return freq_cminverse[0:spectra.size//2], spectra[0:spectra.size//2]

    def fft3(self, x ):
        # Adding zeros to the end of x
        lineshape = fftpack.dct(x, type=1)
        freq_au = np.linspace(0, 0.5/self.dtfs * 1e15, len(x))
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
        if os.path.isfile(local_filename):
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

    def output_Ovelocity_autocorrelation(self):
        local_filename = "%s.Ovac.txt" %self.xyz_filename
        if os.path.isfile(local_filename):
            print("Have calculated Oxygen diffusion for %s, skipping..." %self.xyz_filename)
        else:
            self.cacl_Ovelocity_water_traj()
            # output data
            data = np.zeros((np.size(self.vacf_time_fs_O), 10))
            data[:, 0] = self.vacf_time_fs_O
            data[:, 1] = self.vacf_Ox
            data[:, 2] = self.vacf_Oy
            data[:, 3] = self.vacf_Oz
            data[:, 4] = self.vacf_Otot
            data[:, 5] = self.vacf_Ox_freq
            data[:, 6] = smooth(self.vacf_Ox_sp)
            data[:, 7] = smooth(self.vacf_Oy_sp)
            data[:, 8] = smooth(self.vacf_Oz_sp)
            data[:, 9] = smooth(self.vacf_Otot_sp)
            comments = "# vacf_time_fs, vacf_x, vacf_y, vacf_z, vacf_tot, freq, sp_x, sp_y, sp_z, sp_tot"
            np.savetxt(local_filename, data, comments=comments)

    def output_Hvelocity_autocorrelation(self):
        local_filename = "%s.Hvac.txt" %self.xyz_filename
        if os.path.isfile(local_filename):
            print("Have calculated Hydrogen diffusion for %s, skipping..." %self.xyz_filename)
        else:
            self.cacl_Hvelocity_water_traj()
            # output data
            data = np.zeros((np.size(self.vacf_time_fs_H), 10))
            data[:, 0] = self.vacf_time_fs_H
            data[:, 1] = self.vacf_Hx
            data[:, 2] = self.vacf_Hy
            data[:, 3] = self.vacf_Hz
            data[:, 4] = self.vacf_Htot
            data[:, 5] = self.vacf_Hx_freq
            data[:, 6] = smooth(self.vacf_Hx_sp)
            data[:, 7] = smooth(self.vacf_Hy_sp)
            data[:, 8] = smooth(self.vacf_Hz_sp)
            data[:, 9] = smooth(self.vacf_Htot_sp)
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

    def output_OHkinetic_autocorrelation(self):
        local_filename = "%s.OHkac.txt" %self.xyz_filename
        if os.path.isfile(local_filename):
            print("Have calculated OH bond kinetic for %s, skipping..." %self.xyz_filename)
        else:
            self.cacl_OHkinetic_water_traj()
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

    def output_dipole_autocorrelation(self):
        local_filename = "%s.dac.txt" %self.xyz_filename
        if os.path.isfile(local_filename):
            print("Have calculated dipole autocorrelation for %s, skipping..." %self.xyz_filename)
        else:
            self.cacl_dipole_water_traj()
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
            comments = "# dacf_time_fs, dacf_x, dacf_y, dacf_z, dacf_tot, freq, sp_x, sp_y, sp_z, sp_tot"
            np.savetxt(local_filename, data, comments=comments)

    def output_orientation_autocorrelation(self):
        local_filename = "%s.oac1.txt" %self.xyz_filename
        if os.path.isfile(local_filename):
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

    def output_orientation_autocorrelation2(self):
        local_filename = "%s.oac2.txt" %self.xyz_filename
        if os.path.isfile(local_filename):
            print("Have calculated 2nd orientation autocorrelation for %s, skipping..." %self.xyz_filename)
        else:
            self.cacl_orientation_water_traj2()
            # output data
            data = np.zeros((np.size(self.oacfz2_time_fs), 4))
            data[:, 0] = self.oacfz2_time_fs
            data[:, 1] = self.oacfz2_tot
            data[:, 2] = self.oacfz2_tot_freq
            data[:, 3] = smooth(self.oacfz2_tot_sp)
            comments = "# oacfz2_time_fs, oacfz2_tot, freq,  sp_tot"
            np.savetxt(local_filename, data, comments=comments)

    def output_bond_length(self, bin_start=0.5, bin_end=1.5, bin_num=1000):
        local_filename = "%s.bond_length_dist.txt" %self.xyz_filename
        if os.path.isfile(local_filename):
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
        if os.path.isfile(local_filename):
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

    def output_tot(self):
        self.output_bond_length()
        self.output_pair_dist()
        self.output_velocity_autocorrelation()
        self.output_dipole_autocorrelation()
        self.output_orientation_autocorrelation()
        self.output_orientation_autocorrelation2()



if __name__ == "__main__":
    default_file_size=313681365
    for path in sys.argv[1:]:
        filenames = glob.glob("%s/simu_*.xc.xyz" %path)
        for filename in filenames:
            if os.path.getsize(filename) == default_file_size:
                #a = MD_Analysis(xyz_filename=filename)
                #a.output_bond_length()
                #a = MD_Analysis(xyz_filename=filename)
                #a.output_pair_dist()
                #a = MD_Analysis(xyz_filename=filename)
                #a.output_velocity_autocorrelation()
                a = MD_Analysis(xyz_filename=filename)
                a.output_dipole_autocorrelation()
                #a = MD_Analysis(xyz_filename=filename)
                #a.output_orientation_autocorrelation()
                #a = MD_Analysis(xyz_filename=filename)
                #a.output_orientation_autocorrelation2()
                #a = MD_Analysis(xyz_filename=filename)
                #a.output_Ovelocity_autocorrelation()
                #a = MD_Analysis(xyz_filename=filename)
                #a.output_OHvelocity_autocorrelation()
                #a = MD_Analysis(xyz_filename=filename)
                #a.output_Hvelocity_autocorrelation()
                #a = MD_Analysis(xyz_filename=filename)
                #a.output_OHkinetic_autocorrelation()
