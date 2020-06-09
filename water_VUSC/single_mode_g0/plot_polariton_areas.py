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

def UP_LP_func(x, k):
    # the first half denotes the LP, last half denotes the UP
    size = np.size(x) // 2
    x1 = x[0:size]
    x2 = x[size:]
    omega_0 = 3550
    omega_c = 3550
    OmegaN = k * x1
    LP = np.sqrt((omega_0**2 + omega_c**2 + OmegaN**2)/2.0 - np.sqrt((omega_c**2 + omega_0**2 + OmegaN**2)**2 - 4.0 * omega_0**2 * omega_c**2) / 2.0)
    UP = np.sqrt((omega_0**2 + omega_c**2 + OmegaN**2)/2.0 + np.sqrt((omega_c**2 + omega_0**2 + OmegaN**2)**2 - 4.0 * omega_0**2 * omega_c**2) / 2.0)
    return np.concatenate((LP, UP))

def UP_LP_func_OmegaN(x, OmegaN):
    # the first half denotes the LP, last half denotes the UP
    size = np.size(x) // 2
    x1 = x[0:size]
    x2 = x[size:]
    omega_0 = 3550
    omega_c = 3550
    LP = np.sqrt((omega_0**2 + omega_c**2 + OmegaN**2)/2.0 - np.sqrt((omega_c**2 + omega_0**2 + OmegaN**2)**2 - 4.0 * omega_0**2 * omega_c**2) / 2.0)
    UP = np.sqrt((omega_0**2 + omega_c**2 + OmegaN**2)/2.0 + np.sqrt((omega_c**2 + omega_0**2 + OmegaN**2)**2 - 4.0 * omega_0**2 * omega_c**2) / 2.0)
    return np.concatenate((LP, UP))

def UP_LP_area_OmegaN(x, OmegaN):
    # the first half denotes the LP, last half denotes the UP
    size = np.size(x) // 2
    x1 = x[0:size]
    x2 = x[size:]
    omega_0 = 3550
    omega_c = 3550
    LP = np.sqrt((omega_0**2 + omega_c**2 + OmegaN**2)/2.0 - np.sqrt((omega_c**2 + omega_0**2 + OmegaN**2)**2 - 4.0 * omega_0**2 * omega_c**2) / 2.0)
    UP = np.sqrt((omega_0**2 + omega_c**2 + OmegaN**2)/2.0 + np.sqrt((omega_c**2 + omega_0**2 + OmegaN**2)**2 - 4.0 * omega_0**2 * omega_c**2) / 2.0)
    theta = np.arctan2(2.0*omega_c*OmegaN, OmegaN**2 + omega_0**2 - omega_c**2)
    theta[0] = np.pi/2.0
    print("OmegaN is", OmegaN)
    print("theta is", theta)
    return np.concatenate((LP**2 * np.sin(theta/2.0)**2, UP**2 * np.cos(theta/2.0)**2))


from scipy.optimize import curve_fit

def plot_polaritons():
    paths = ["E0_0e-4", "E0_1e-4", "E0_2e-4", "E0_3e-4", "E0_4e-4", "E0_5e-4", "E0_6e-4", "E0_7e-4", "E0_8e-4"]
    LPs = np.array([3550,         3481,      3417,       3351,       3289,    3249,      3212,       3181,       3153])
    UPs = np.array([3550,         3685,      3817,       3916,       4030,    4177,      4334,       4492,       4662])
    x0 = np.array([0.0, 1e-4, 2e-4, 3e-4, 4e-4, 5e-4, 6e-4, 7e-4, 8e-4])
    LPs_area, UPs_area, LPs_freq, UPs_freq = [], [], [], []
    for i, path in enumerate(paths):
        data = obtain_avg_data(path, pattern="simu_*.dac.txt")
        x, y = data[:,5], (data[:,6] + data[:,7])/2e28
        dx = x[2] - x[1]
        # locate the index of LP and UP
        index1 = int(LPs[i] / dx)
        index2 = int(UPs[i] / dx)
        if i is 0:
            minus, plus, minus2, plus2 = 250, 250, 250, 250
        elif i is 1:
            minus, plus, minus2, plus2 = 100, 59, 61, 150
        else:
            minus, plus, minus2, plus2 = 130, 100, 130, 100
        LP = y[index1-minus:index1+plus]
        LP_freq = x[index1-minus:index1+plus]
        UP = y[index2-minus2:index2+plus2]
        UP_freq = x[index2-minus2:index2+plus2]
        LParea = np.sum(LP)*dx
        UParea = np.sum(UP)*dx
        LPs_area.append(LParea if i > 0 else LParea/2.0)
        UPs_area.append(UParea if i > 0 else UParea/2.0)
        LPs_freq.append(LP_freq[np.argmax(LP)])
        UPs_freq.append(UP_freq[np.argmax(UP)])
    LPs_area = np.array(LPs_area)
    UPs_area = np.array(UPs_area)
    LPs_freq = np.array(LPs_freq)
    UPs_freq = np.array(UPs_freq)

    axes = clp.initialize(3, 1, width=4.3, height=9.5, LaTeX=True, fontsize=12,  labelthem=True, labelthemPosition=[0.1, 0.95])
    xs = [x0]*2
    Rabi_freq = UPs_freq - LPs_freq
    xs2 = [Rabi_freq / 2.0 / 3550.0] * 2
    ys1 = [LPs_freq, UPs_freq]
    ys2 = [LPs_area, UPs_area]
    labels = ["LP", "UP"]
    colors = ["b*", "ro"]

    popt, pocv = curve_fit(linear_func, x0, Rabi_freq)
    print("popt is %.4E" %popt[0])

    for i in range(3):
        if i is 0:
            axes[i].axvspan(0, 1.91e-4*2.0, alpha=0.1, color='g')
            axes[i].axvspan(1.91e-4*2.0, 8.3e-4, alpha=0.1, color='r')
        axes[i].text(0.2, 0.9, "VSC", transform=axes[i].transAxes, fontsize=12, color="g")
        axes[i].text(0.54, 0.9, "V-USC", transform=axes[i].transAxes, fontsize=12, color="r")

    for i in range(2):
        axes[i+1].axvspan(0, 0.1, alpha=0.1, color='g')
        axes[i+1].axvspan(0.1, 8.3e-4*popt[0], alpha=0.1, color='r')

    clp.plotone([x0]*2, [Rabi_freq, linear_func(x0, *popt)], axes[0], colors=["k^", "0.5"], labels=["simulation", "linear fit"],
            ylabel="Rabi frequency [cm$^{-1}$]", xlim=[0, 8.3e-4], xlabel="$\widetilde{\\varepsilon}$ [a.u.]")
    xmax = 8.3e-4 * popt[0] / 2.0 / 3550.0
    print(xmax)


    clp.plotone(xs2, ys1, axes[1], colors=colors, labels=labels, ylabel="Polariton frequency [cm$^{-1}$]", xlim=[0, xmax])
    x_target = np.concatenate((x0, x0))
    y_target = np.concatenate((LPs_freq, UPs_freq))
    #popt, pocv = curve_fit(UP_LP_func, x_target, y_target, p0=[2e6])
    results = UP_LP_func_OmegaN(x_target, Rabi_freq)
    size = np.size(results)//2
    clp.plotone([xs2[0]], [results[0:size]], axes[1], colors=["k--"], labels=["1D analytical"], xlim=[0, xmax],xlabel="$\Omega_N / 2\omega_0$")
    clp.plotone([xs2[0]], [results[size:]], axes[1], colors=["k--"], showlegend=False, xlim=[0, xmax])

    clp.plotone(xs2, ys2, axes[2], colors=colors, labels=labels, ylabel="Integrated peak area [arb. units]",
                xlabel="$\Omega_N / 2\omega_0$", xlim=[0, xmax], showlegend=False)
    results = UP_LP_area_OmegaN(x_target, Rabi_freq)
    LP, UP = results[0:size], results[size:]
    clp.plotone([xs2[0]], [ LP / LP[0] * LPs_area[0]], axes[2], colors=["k--"], xlim=[0, xmax], showlegend=False)
    clp.plotone([xs2[0]], [ UP / UP[0] * UPs_area[0]], axes[2], colors=["k--"], showlegend=False, xlim=[0, xmax])

    axes[0].legend(loc="center right")
    axes[1].legend(loc="center right")
    #axes[2].legend(loc="center right")
    clp.adjust(savefile="Rabi_splitting.pdf")

if __name__ == "__main__":
    plot_polaritons()
