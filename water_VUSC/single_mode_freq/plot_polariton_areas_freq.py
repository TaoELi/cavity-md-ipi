import numpy as np
import columnplots as clp
from scipy import signal
import os, sys, time
import math
from itertools import islice
from scipy import fftpack
import glob
import json

def obtain_avg_data(path, pattern="simu_*.vac.txt"):
    filenames = glob.glob("%s/%s" %(path, pattern))
    print("Averaging over %d trajectories" %(len(filenames)))
    data = np.loadtxt(filenames[0])
    for filename in filenames[1:]:
        data += np.loadtxt(filename)
    data /= float(len(filenames))
    return data

def linear_func(x, k):
    return k*x

from scipy.optimize import curve_fit

def calc_OmegaN(omega_c):
    # the first half denotes the LP, last half denotes the UP
    omega_0 = 3550
    OmegaN = 927
    LP = np.sqrt((omega_0**2 + omega_c**2 + OmegaN**2)/2.0 - np.sqrt((omega_c**2 + omega_0**2 + OmegaN**2)**2 - 4.0 * omega_0**2 * omega_c**2) / 2.0)
    UP = np.sqrt((omega_0**2 + omega_c**2 + OmegaN**2)/2.0 + np.sqrt((omega_c**2 + omega_0**2 + OmegaN**2)**2 - 4.0 * omega_0**2 * omega_c**2) / 2.0)
    return UP - LP

def calc_LP_UP(omega_c):
    # the first half denotes the LP, last half denotes the UP
    omega_0 = 3550
    OmegaN = 927
    LP = np.sqrt((omega_0**2 + omega_c**2 + OmegaN**2)/2.0 - np.sqrt((omega_c**2 + omega_0**2 + OmegaN**2)**2 - 4.0 * omega_0**2 * omega_c**2) / 2.0)
    UP = np.sqrt((omega_0**2 + omega_c**2 + OmegaN**2)/2.0 + np.sqrt((omega_c**2 + omega_0**2 + OmegaN**2)**2 - 4.0 * omega_0**2 * omega_c**2) / 2.0)
    return LP, UP

def calc_LP_UP_area(omega_c):
    # the first half denotes the LP, last half denotes the UP
    omega_0 = 3550
    OmegaN = 927
    LP = np.sqrt((omega_0**2 + omega_c**2 + OmegaN**2)/2.0 - np.sqrt((omega_c**2 + omega_0**2 + OmegaN**2)**2 - 4.0 * omega_0**2 * omega_c**2) / 2.0)
    UP = np.sqrt((omega_0**2 + omega_c**2 + OmegaN**2)/2.0 + np.sqrt((omega_c**2 + omega_0**2 + OmegaN**2)**2 - 4.0 * omega_0**2 * omega_c**2) / 2.0)
    theta = np.arctan2(2.0*omega_c*OmegaN, OmegaN**2 + omega_0**2 - omega_c**2)
    return LP**2 * np.sin(theta/2.0)**2, UP**2 * np.cos(theta/2.0)**2

def find_cross_point(OmegaN):
    omega_0 = 3550
    omega_c = np.linspace(3000, 6000, 1000)
    LP = np.sqrt((omega_0**2 + omega_c**2 + OmegaN**2)/2.0 - np.sqrt((omega_c**2 + omega_0**2 + OmegaN**2)**2 - 4.0 * omega_0**2 * omega_c**2) / 2.0)
    UP = np.sqrt((omega_0**2 + omega_c**2 + OmegaN**2)/2.0 + np.sqrt((omega_c**2 + omega_0**2 + OmegaN**2)**2 - 4.0 * omega_0**2 * omega_c**2) / 2.0)
    theta = np.arctan2(2.0*omega_c*OmegaN, OmegaN**2 + omega_0**2 - omega_c**2)
    result = np.abs(LP**2 * np.sin(theta/2.0)**2 - UP**2 * np.cos(theta/2.0)**2)
    return omega_c[np.argmin(result)]

def calc_cross_points():
    OmegaN = np.linspace(0.0, 2000, 100)
    return OmegaN, np.array([find_cross_point(x) for x in OmegaN])

def plot_polaritons():
    paths = ["Freq_%d" %(3150+200*i) for i in range(12)]
    LPs = [3017, 3174, 3258, 3353, 3416, 3437, 3468, 3489, 3510, 3500, 3521, 3521]
    UPs = [4088, 4120, 4183, 4288, 4393, 4519, 4666, 4823, 4991, 5180, 5348, 5537]
    x0 = np.array([3150+200*i for i in range(12)])
    LPs_area, UPs_area, LPs_freq, UPs_freq = [], [], [], []
    for i, path in enumerate(paths):
        data = obtain_avg_data(path, pattern="simu_*.dac.txt")
        x, y = data[:,5], (data[:,6] + data[:,7])/2e28
        dx = x[2] - x[1]
        # locate the index of LP and UP
        index1 = int(LPs[i] / dx)
        index2 = int(UPs[i] / dx)
        minus, plus, minus2, plus2 = 150, 300, 150, 150
        LP = y[index1-minus:index1+plus]
        LP_freq = x[index1-minus:index1+plus]
        UP = y[index2-minus2:index2+plus2]
        UP_freq = x[index2-minus2:index2+plus2]
        LParea = np.sum(LP)*dx
        UParea = np.sum(UP)*dx
        LPs_area.append(LParea)
        UPs_area.append(UParea)
        LPs_freq.append(LP_freq[np.argmax(LP)])
        UPs_freq.append(UP_freq[np.argmax(UP)])
    LPs_area = np.array(LPs_area)
    UPs_area = np.array(UPs_area)
    LPs_freq = np.array(LPs_freq)
    UPs_freq = np.array(UPs_freq)

    axes = clp.initialize(2, 1, width=4.3, height=6.0, LaTeX=True, fontsize=12, sharex=True, labelthem=True, labelthemPosition=[0.1, 0.95])
    xs = [x0]*2
    ys1 = [LPs_freq, UPs_freq]
    ys2 = [LPs_area, UPs_area]
    labels = ["LP", "UP"]
    colors = ["b*", "ro"]
    clp.plotone(xs, ys1, axes[0], colors=colors, labels=labels, ylabel="Polariton frequency [cm$^{-1}$]")
    LP_calc, UP_calc = calc_LP_UP(x0)
    clp.plotone([x0], [LP_calc], axes[0], colors=["k--"], labels=["1d analytical"])
    clp.plotone([x0], [UP_calc], axes[0], colors=["k--"], showlegend=False)
    clp.plotone([x0]*2, [3550 * np.ones(np.shape(x0)), x0], axes[0], colors=["0.5", "0.5"], showlegend=False)

    #clp.plotone([x0], [3550 * np.ones(np.shape(x0))], axes[0], colors="k--", showlegend=False)
    clp.plotone(xs, ys2, axes[1], colors=colors,  ylabel="Integrated peak area [arb. units]", xlabel="cavity mode frequency [cm$^{-1}$]", showlegend=False, ylim=[-3, 370])
    LP_calc, UP_calc = calc_LP_UP_area(x0)
    clp.plotone([x0], [LP_calc * LPs_area[5] / LP_calc[5]], axes[1], colors=["k--"], showlegend=False)
    clp.plotone([x0], [UP_calc * UPs_area[5] / UP_calc[5]], axes[1], colors=["k--"], showlegend=False)

    from mpl_toolkits.axes_grid.inset_locator import inset_axes
    xsmall, ysmall = calc_cross_points()
    inset_ax1 = inset_axes(axes[1], width="40%", height=0.6, loc=1)
    clp.plotone([xsmall / 3550.0 / 2.0], [ysmall], inset_ax1, colors=["k--"], showlegend=False, lw=1.50)#, xlabel="Omega_N / 2\omega_0", ylabel="cross point frequency [cm$^{-1}$]")
    axes[1].text(4700, 198, "$\Omega_N / 2\omega_0$")
    axes[1].text(4000, 180, "cross freq. [cm$^{-1}$]", rotation=90)

    axes[0].legend(loc="center right")
    #axes[2].legend(loc="center right")
    clp.adjust(savefile="Rabi_splitting_freq.pdf")

if __name__ == "__main__":
    plot_polaritons()
