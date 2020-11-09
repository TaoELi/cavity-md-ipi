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
import matplotlib.patches as patches


def obtain_avg_data(path, pattern="simu_*.vac.txt"):
    filenames = glob.glob("%s/%s" %(path, pattern))
    print("Averaging over %d trajectories under %s" %(len(filenames), path))
    data = np.loadtxt(filenames[0])
    for filename in filenames[1:]:
        data += np.loadtxt(filename)
    data /= float(len(filenames))
    return data

def prepare_data(path="exciteCO2UP", pattern="simu_*.localIR_td.txt", freq_start=2100, freq_end=2500):
    data2 = obtain_avg_data(path=path+"/Amp_6e-3", pattern=pattern)
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


from scipy.optimize import curve_fit
def fit_exponential(x, y):
    popt, pocv =  curve_fit(func, x, y)
    yfit = func(x, *popt)
    print("popt is", popt)
    return yfit, popt[0]

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

def obtain_polariton_lifetime(path, omega=2320, nstart=300):
    t, phx, phy = prepare_ph_data(path=path, omega=omega)
    yfit_phx, tau = fit_exponential(t[nstart:], phx[nstart:] + phy[nstart:])
    return tau, t, phx + phy

def plot_IR():
    freq, data2_LP_tot = prepare_data(path="exciteCO2LP", pattern="simu_*.dacf_td.txt")
    freq, data2_LP = prepare_data(path="exciteCO2LP", pattern="simu_*.localIR_td.txt")
    freq, data2_LP_tot_nocav = prepare_data(path="exciteCO2LP_nocavity", pattern="simu_*.dacf_td.txt")
    freq, data2_LP_nocav = prepare_data(path="exciteCO2LP_nocavity", pattern="simu_*.localIR_td.txt")
    freq, data2_LP_tot_2200 = prepare_data(path="exciteCO2LP_2200", pattern="simu_*.dacf_td.txt")
    freq, data2_LP_2200 = prepare_data(path="exciteCO2LP_2200", pattern="simu_*.localIR_td.txt")

    ratio = 4.0

    extent = [0, 60, freq[0] , freq[-1]]
    from matplotlib.colors import LogNorm, Normalize
    vmax1 = np.max(np.max(data2_LP_tot)) * 1.0
    vmin1 = vmax1 * 3e-3
    vmax2 = np.max(np.max(data2_LP)) * 1.0
    vmin2 = vmax2 * 0.0
    fig, axes = clp.initialize(2, 3, width=4.3*2.2, height=4.0, LaTeX=True, fontsize=12, return_fig_args=True, labelthem=True, labelthemPosition=[-0.03, 1.12])
    im1 = axes[0, 0].imshow(data2_LP_tot, aspect='auto', extent=extent,  cmap=cm.inferno, interpolation=None, norm=LogNorm(vmin=vmin1, vmax=vmax1))
    im1 = axes[0, 1].imshow(data2_LP_tot_2200, aspect='auto', extent=extent,  cmap=cm.inferno, interpolation=None, norm=LogNorm(vmin=vmin1, vmax=vmax1))
    im1 = axes[0, 2].imshow(data2_LP_tot_nocav * ratio, aspect='auto', extent=extent,  cmap=cm.inferno, interpolation=None, norm=LogNorm(vmin=vmin1, vmax=vmax1))
    im2 = axes[1, 0].imshow(data2_LP, aspect='auto', extent=extent,  cmap=cm.hot, interpolation=None, norm=Normalize(vmin=vmin2, vmax=vmax2))
    im2 = axes[1, 1].imshow(data2_LP_2200, aspect='auto', extent=extent,  cmap=cm.hot, interpolation=None, norm=Normalize(vmin=vmin2, vmax=vmax2))
    im2 = axes[1, 2].imshow(data2_LP_nocav * ratio, aspect='auto', extent=extent,  cmap=cm.hot, interpolation=None, norm=Normalize(vmin=vmin2, vmax=vmax2))

    x = range(0,60)
    xs = [x]*2
    ys = [np.ones(len(x)) * 2320, np.ones(len(x)) * 2260]
    clp.plotone(xs, ys, axes[0,0], showlegend=False, colors=["w--", "c--"], lw=1., ylabel="IR freq [cm$^{-1}$]", alpha=0.0)
    clp.plotone(xs, ys, axes[1,0], showlegend=False, colors=["w--", "c--"], lw=1., ylabel="local IR freq [cm$^{-1}$]", xlabel="Time [ps]", alpha=0.0)
    clp.plotone(xs, ys, axes[0,1], showlegend=False, colors=["w--", "c--"], lw=1., alpha=0.0)
    clp.plotone(xs, ys, axes[1,1], showlegend=False, colors=["w--", "c--"], lw=1., xlabel="time [ps]", alpha=0.0)
    clp.plotone(xs, ys, axes[0,2], showlegend=False, colors=["w--", "c--"], lw=1., alpha=0.0)
    clp.plotone(xs, ys, axes[1,2], showlegend=False, colors=["w--", "c--"], lw=1., xlabel="time [ps]", alpha=0.0)

    axes[0, 1].set_yticklabels([])
    axes[1, 1].set_yticklabels([])
    axes[0, 2].set_yticklabels([])
    axes[1, 2].set_yticklabels([])
    axes[0, 0].set_xticklabels([])
    axes[0, 1].set_xticklabels([])
    axes[0, 2].set_xticklabels([])

    axes[0, 0].text(20, 2440, "cavity mode 2320 cm$^{-1}$", color='w')
    axes[0, 1].text(20, 2440, "cavity mode 2200 cm$^{-1}$", color='w')
    axes[0, 2].text(35, 2440, "outside cavity", color='w')
    axes[1, 2].text(5, 2440, "$\\times$ %.1f" %ratio, color='w', fontsize=16)
    axes[0, 2].text(5, 2440, "$\\times$ %.1f" %ratio, color='w', fontsize=16)

    arrow_properties = dict(facecolor="r", width=0.5, headwidth=4, shrink=0.1)
    pulse_locations = [2241, 2167, 2241]
    for i in range(2):
        for j in range(3):
            axes[i, j].tick_params(color='c', labelsize='medium', width=2)
            axes[i, j].annotate("", xy=(0, pulse_locations[j]), xytext=(-6, pulse_locations[j]), arrowprops=arrow_properties)

    fig.subplots_adjust(right=0.95)
    cbar_ax = fig.add_axes([0.96, 0.53, 0.02, 0.35])
    fig.colorbar(im1, cax=cbar_ax)
    cbar_ax = fig.add_axes([0.96, 0.11, 0.02, 0.35])
    fig.colorbar(im2, cax=cbar_ax)

    for i in range(3):
        rect = patches.Rectangle((0,2150), 60, 70, linewidth=1., edgecolor='c',facecolor='none')
        axes[1, i].add_patch(rect)
        rect = patches.Rectangle((0,2225), 60, 140, linewidth=1.3, edgecolor='0.5',facecolor='none')
        axes[1, i].add_patch(rect)
    
    fig.suptitle('Time-resolved spectra after exciting LP', fontsize=14, weight='extra bold')

    clp.adjust(savefile="dynamics_CO2_LP.pdf")

    # Finally, we fit the exponential to calculate the lifetime of the polaritons
    freq, data2_LP = prepare_data(path="exciteCO2LP", pattern="simu_*.dacf_td.txt", freq_start=2200, freq_end=2280)
    y = np.sum(data2_LP, axis=0)
    N = np.size(y)
    x = np.linspace(0, 60, N)
    from scipy.optimize import curve_fit
    def func(x, k, a, b):
        return a * np.exp(-x/k) + b
    popt, pcov = curve_fit(func, x, y)
    yfit = func(x, *popt)
    print("popt is", popt)
    ax = clp.initialize()
    clp.plotone([x]*2, [y, yfit], ax, colors=["ko", "r--"], labels=["LP total intensities", "fit lifetime %.1f ps" %popt[0]])
    clp.adjust(savefile="summed_excitation_statistics_CO2_LP.pdf")

    # Second, we fit the exponential to calculate the lifetime of the polaritons
    freq, data2_LP = prepare_data(path="exciteCO2LP", pattern="simu_*.dacf_td.txt", freq_start=2380, freq_end=2480)
    y = np.sum(data2_LP, axis=0)
    N = np.size(y)
    x = np.linspace(0, 60, N)
    from scipy.optimize import curve_fit
    def func(x, k, a, b):
        return a * np.exp(-x/k) + b
    popt, pcov = curve_fit(func, x, y)
    yfit = func(x, *popt)
    print("popt is", popt)
    ax = clp.initialize()
    clp.plotone([x]*2, [y, yfit], ax, colors=["ko", "r--"], labels=["UP total intensities", "fit lifetime %.1f ps" %popt[0]])
    clp.adjust(savefile="summed_excitation_statistics_CO2_UP.pdf")

    freq, data2_LP = prepare_data(path="exciteCO2LP_2200", pattern="simu_*.dacf_td.txt", freq_start=2100, freq_end=2200)
    y = np.sum(data2_LP, axis=0)
    N = np.size(y)
    x = np.linspace(0, 60, N)
    from scipy.optimize import curve_fit
    def func(x, k, a, b):
        return a * np.exp(-x/k) + b
    popt, pcov = curve_fit(func, x, y)
    yfit = func(x, *popt)
    print("popt is", popt)
    ax = clp.initialize()
    clp.plotone([x]*2, [y, yfit], ax, colors=["ko", "r--"], labels=["LP total intensities", "fit lifetime %.1f ps" %popt[0]])
    clp.adjust(savefile="summed_excitation_statistics_CO2_LP_2200.pdf")

    # Finally, we fit the exponential to calculate the lifetime of the polaritons
    freq, data2_LP = prepare_data(path="exciteCO2LP_2200", pattern="simu_*.dacf_td.txt", freq_start=2350, freq_end=2450)
    y = np.sum(data2_LP, axis=0)
    N = np.size(y)
    x = np.linspace(0, 60, N)
    from scipy.optimize import curve_fit
    def func(x, k, a, b):
        return a * np.exp(-x/k) + b
    popt, pcov = curve_fit(func, x, y)
    yfit = func(x, *popt)
    print("popt is", popt)
    ax = clp.initialize()
    clp.plotone([x]*2, [y, yfit], ax, colors=["ko", "r--"], labels=["UP total intensities", "fit lifetime %.1f ps" %popt[0]])
    clp.adjust(savefile="summed_excitation_statistics_CO2_UP_2200.pdf")

    # also fit the decay rate of LP under 2200 with photonic dynamics
    tau, t, ph = obtain_polariton_lifetime(path="exciteCO2LP_2200/Amp_6e-3", omega=2200, nstart=300)
    ax = clp.initialize()
    clp.plotone([t], [ph], ax, colors=["ko"], labels=["ph dynamics lifetime %.1f ps" %tau])
    clp.adjust(savefile="summed_excitation_statistics_CO2_LP_2200_ph.pdf")

if __name__ == "__main__":
    plot_IR()
