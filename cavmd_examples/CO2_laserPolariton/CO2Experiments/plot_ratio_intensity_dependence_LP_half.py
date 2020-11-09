import numpy as np
import columnplots as clp
from scipy import signal
import os, sys, time
import math
from itertools import islice
from scipy import fftpack
import glob
import json
import matplotlib.cm as cm
import matplotlib.transforms as mtransforms

UseUnitAU=False

def E0toE0(x, UseUnitAU=UseUnitAU):
    if UseUnitAU:
        return x
    else:
        return x**2 / 6e-3 /6e-3 * 632

tmax = 60

def obtain_avg_data(path, pattern="simu_*.vac.txt"):
    filenames = glob.glob("%s/%s" %(path, pattern))
    print("Averaging over %d trajectories under %s" %(len(filenames), path))
    data = np.loadtxt(filenames[0])
    for filename in filenames[1:]:
        data += np.loadtxt(filename)
    data /= float(len(filenames))
    return data

def prepare_ph_data(path="exciteCO2LP/Amp_6e-3", omega=2320, pattern="simu_*ph*.txt"):
    filenames = glob.glob("%s/%s" %(path, pattern))
    N = len(filenames)
    print("reading %s with %d files" %(path, N))
    size = np.size(np.loadtxt(filenames[0])[:,0])
    t, phx, phy = np.zeros(size), np.zeros(size), np.zeros(size)
    t2ps = 1e-3
    ph2au = 0.5 * 1.8897259885789**2 * omega**2 * 0.0000045563352812122295**2
    vph2au = 0.5 * 1.8897259885789**2 / 41.341374575751**2
    for filename in filenames:
        data = np.loadtxt(filename)
        t += data[:,0]*t2ps
        phx += data[:,1]**2 * ph2au + data[:,3]**2 * vph2au
        phy += data[:,2]**2 * ph2au + data[:,4]**2 * vph2au
    return t/N, phx/N, phy/N

def prepare_data(path="exciteCO2UP", pattern="simu_*.localIR_td.txt", freq_start=2100, freq_end=2500):
    data2 = obtain_avg_data(path=path, pattern=pattern)
    freq, data2 = data2[:, 0], data2[:, 1:60]
    x = freq
    dx = x[2] - x[1]
    nstart, nend = int(freq_start / dx), int(freq_end / dx)
    x = x[nstart:nend]
    data2 = data2[nstart:nend,:]
    data2 = np.abs((data2[::-1, :]) / 2e28)
    return x, data2

def func(x, k, a, b):
    return a * np.exp(-x/k) + b

def func2(x, k, k1, a, b, c, t0):
    return (a * np.exp(-(x-t0)/k) + b * np.exp(-(x-t0)/k1) + c) * (x > t0) + (x <= t0) * (a + b + c)

def func_linear(x, k, b):
    return k*x + b

def func_quadratic(x, k2, b):
    return k2*x**2 + b

from scipy.optimize import curve_fit
def fit_exponential(x, y):
    popt, pocv =  curve_fit(func, x, y)
    yfit = func(x, *popt)
    print("popt is", popt)
    return yfit, popt[0]

def fit_biexponential(x, y):
    popt, pocv =  curve_fit(func2, x, y, maxfev=12000)
    yfit = func2(x, *popt)
    print("popt is", popt)
    return yfit, popt[0], popt[1]

def fit_linear(x, y):
    popt, pocv =  curve_fit(func_linear, x, y, maxfev=12000)
    yfit = func_linear(x, *popt)
    print("popt is", popt)
    return yfit, popt[0]

def fit_quadratic(x, y):
    popt, pocv =  curve_fit(func_quadratic, x, y, maxfev=12000)
    yfit = func_quadratic(x, *popt)
    print("popt is", popt)
    return yfit, popt

def plot_photon():
    amps = ['6e-4', '8e-4', '1e-3', '1.5e-3', '2e-3', '3e-3', '4e-3', '5e-3', '6e-3', '7e-3', '8e-3', '9e-3']
    paths = ["exciteCO2LP/Amp_%s" %amp for amp in amps]
    paths_nocav = ["exciteCO2LP_nocavity/Amp_%s" %amp for amp in amps]
    print("Paths are", paths)
    amps = [float(amp) for amp in amps]
    x = np.array(amps, dtype=np.float32)

    # Prepare data for photonic decay lifetime as a function of amplitudes
    '''
    y_ph = []
    for path in paths:
        t, phx, phy = prepare_ph_data(path=path)
        try:
            yfit_phx, tau = fit_exponential(t[300:], phx[300:]+phy[300:])
            print(tau)
            y_ph.append(tau)
        except:
            None
    y_ph = np.array(y_ph)
    '''

    # Prepare data for v=2/v=1 intensities as a function of amplitudes
    y_ratio = []
    for path in paths:
        #freq, data2_LP_v1 = prepare_data(path=path, freq_start=2220, freq_end=2360)
        #y_v1 = np.sum(data2_LP_v1, axis=0)
        freq, data2_LP_v2 = prepare_data(path=path, freq_start=2150, freq_end=2220)
        y_v2 = np.sum(data2_LP_v2, axis=0)
        y = y_v2  #/ y_v1
        #y_ratio.append(np.max(y))
        y_ratio.append(y[1])
    y_ratio = np.array(y_ratio)

    y_ratio_nocav = []
    for path in paths_nocav:
        #freq, data2_LP_v1 = prepare_data(path=path, freq_start=2220, freq_end=2360)
        #y_v1 = np.sum(data2_LP_v1, axis=0)
        freq, data2_LP_v2 = prepare_data(path=path, freq_start=2150, freq_end=2220)
        y_v2 = np.sum(data2_LP_v2, axis=0)
        y = y_v2 #/ y_v1
        #y_ratio_nocav.append(np.max(y))
        y_ratio_nocav.append(y[1])
    y_ratio_nocav = np.array(y_ratio_nocav)

    ax = clp.initialize(1,1,width=4.3, height=4.3*0.618, LaTeX=True, fontsize=12)

    color = 'b'

    x =  E0toE0(x)

    # Performing fit
    # Fit to obtain the slope
    #yfit_linear, k = fit_linear(np.log(x[-5:]), np.log(y_ratio[-5:]))
    yfit_quadratic_cav, k = fit_quadratic(x[5:], y_ratio[5:])
    print("Coeff is", k)

    #yfit_linear, k = fit_linear(np.log(x[-5:]), np.log(y_ratio_nocav[-5:]))
    yfit_quadratic_nocav, k = fit_quadratic(x[5:], y_ratio_nocav[5:])
    print("Coeff is", k)

    #clp.plotone([x], [1.0/y_ph], axes[0], colors=[color+"-^"], ylabel="LP decay rate [ps$^{-1}$]",  showlegend=False,  xlim=[E0toE0(5e-4), E0toE0(9.5e-3)], xlog=True, ylog=True)

    #clp.plotone([x[5:], x[5:], x, x], [yfit_quadratic_cav, yfit_quadratic_nocav, y_ratio,  y_ratio_nocav], axes[1], colors=["0.5", "0.5", color+"o", "ks"], ylabel='Initial-time $I(v_2)$ [arb. units]', xlabel="$E_0$ [a.u.]" if UseUnitAU else "F [mJ/cm$^2$]", showlegend=False)#, xlog=True, ylog=True)
    clp.plotone([x, x], [y_ratio,  y_ratio_nocav], ax, colors=[color+"-o", "k-s"], ylabel='Early-time $I_{NL}$ [arb. units]', xlabel="$E_0$ [a.u.]" if UseUnitAU else "F [mJ/cm$^2$]", showlegend=False, xlog=True, ylog=True, xlim=[E0toE0(5e-4), E0toE0(9.5e-3)], ylim=[8e-5, 9e-1])


    ax.text(0.2, 0.4, "inside cavity", color=color, fontsize=12, transform = ax.transAxes)
    ax.text(0.52, 0.03, "outside cavity", color="k", fontsize=12, transform = ax.transAxes)

    #axes[0].axvspan(E0toE0(4e-4), E0toE0(1.5e-3), alpha=0.2, color='g')
    #axes[0].axvspan(E0toE0(1.5e-3), E0toE0(2e-2), alpha=0.2, color='r')
    ax.axvspan(E0toE0(4e-4), E0toE0(1.5e-3), alpha=0.2, color='g')
    ax.axvspan(E0toE0(1.5e-3), E0toE0(2e-2), alpha=0.2, color='r')
    #axes[0].text(0.1, 0.9, "linear regime", color="g", fontsize=12,  transform = axes[0].transAxes)
    #axes[0].text(0.52, 0.9, "nonlinear regime", color="r", fontsize=12,  transform = axes[0].transAxes)

    clp.adjust(savefile="ratio_LP_E0dependence_half.pdf", tight_layout=True)

if __name__ == "__main__":
    plot_photon()
