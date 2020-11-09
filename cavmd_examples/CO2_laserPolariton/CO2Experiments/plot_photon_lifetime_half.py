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
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

'''
E0: incoming field amplitude in a.u.
change to
E0=6e-3 <-> I = 1.2636e12 W/cm^2 * 5E-13 s =  6.32E-1 J/cm2 = 632 mJ/cm2
E0=6e-4 <-> I = 1.2636e10 W/cm^2 * 5E-13 s =  6.32E-3 J/cm2 = 6.32 mJ/cm2
'''

UseUnitAU=False
E0string_weak = "$E_0=6\\times 10^{-4}$" if  UseUnitAU else "F = 6.32 mJ/cm$^2$"
E0string_strong = "$E_0=6\\times 10^{-3}$" if  UseUnitAU else "F = 632 mJ/cm$^2$"


tmax = 60

UP = np.array([2350, 2375, 2428, 2485, 2548, 2616])
LP = np.array([2301, 2280, 2242, 2207, 2177, 2150])

def obtain_avg_data(path, pattern="simu_*.vac.txt"):
    filenames = glob.glob("%s/%s" %(path, pattern))
    print("Averaging over %d trajectories under %s" %(len(filenames), path))
    data = np.loadtxt(filenames[0])
    for filename in filenames[1:]:
        data += np.loadtxt(filename)
    data /= float(len(filenames))
    return data

def func(x, k, a, b):
    return a * np.exp(-x/k) + b

def func2(x, k, k1, a, b, c, t0):
    return (a * np.exp(-(x-t0)/k) + b * np.exp(-(x-t0)/k1) + c) * (x > t0) + (x <= t0) * (a + b + c)

def func_Lorentz(x, Gamma, A, omega0):
    return A / np.pi * 0.5 * Gamma / ( (x - omega0)**2 + (0.5*Gamma)**2 )

def func_Lorentz_centered(x, Gamma, A):
    return A / np.pi * 0.5 * Gamma / ( (x - 2327)**2 + (0.5*Gamma)**2 )

def func_biLorentz(x, Gamma, A, Gamma2, A2, shift):
    return A / np.pi * 0.5 * Gamma / ( (x - 2327)**2 + (0.5*Gamma)**2 ) + A2 / np.pi * 0.5 * Gamma2 / ( (x - (2251+shift))**2 + (0.5*Gamma2)**2 )

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

def fit_lorentzian(x, y):
    popt, pocv = curve_fit(func_Lorentz, x, y)
    yfit = func_Lorentz(x, *popt)
    print("popt is", popt)
    return yfit, popt[0]

def fit_lorentzian(x, y):
    popt, pocv = curve_fit(func_Lorentz, x, y)
    yfit = func_Lorentz(x, *popt)
    print("popt is", popt)
    return yfit, popt[0]

def fit_Lorentzian_centered(x, y):
    popt, pocv = curve_fit(func_Lorentz_centered, x, y, maxfev=120000)
    yfit = func_Lorentz_centered(x, *popt)
    print("popt is", popt)
    return yfit, popt[0]

def fit_biLorentzian(x, y):
    popt, pocv = curve_fit(func_biLorentz, x, y, maxfev=120000)
    yfit = func_biLorentz(x, *popt)
    print("popt is", popt)
    return yfit, popt[0]

def obtain_peak_linshape(path="CO2only_changeE0/E0_0e-4", freq_start=2220, freq_end=2420):
    data = obtain_avg_data(path=path, pattern="simu_*dac.txt")
    freq, sp = data[:,5], (data[:, 6] + data[:, 7]) / 2e28
    df = freq[1] - freq[0]
    nstart, nend = int(freq_start / df), int(freq_end / df)
    yfit, k1 = fit_lorentzian(freq[nstart:nend], sp[nstart:nend])
    #ax = clp.initialize()
    #clp.plotone([freq, freq[nstart:nend]], [sp, yfit], ax, colors=["r", "b--"], labels=["simu", "fit"], xlim=[0,3000])
    #clp.adjust()
    return k1

def collect_polariton_lineshape():
    couplings = ["5e-5", "1e-4", "2e-4", "3e-4", "4e-4", "5e-4"]
    UP = np.array([2350, 2375, 2428, 2485, 2548, 2616])
    LP = np.array([2301, 2280, 2242, 2207, 2177, 2150])
    # Check if the data is calculated before
    filename = "polariton_lineshape.tmp"
    try:
        data = np.loadtxt(filename)
        print("Have calculated lineshape :)")
        k_bare, kLPs, kUPs = data[1,1], data[:,2], data[:,3]
    except:
        print("From stratch fitting linshapes...")
        # outside cavity
        path ="CO2only_changeE0/E0_0e-4"
        k_bare = obtain_peak_linshape(path=path, freq_start=2120, freq_end=2520)
        kLPs, kUPs = [], []
        for i, c in enumerate(couplings):
            path ="CO2only_changeE0/E0_%s" %c
            k_LP = obtain_peak_linshape(path=path, freq_start=LP[i]-20*i-30, freq_end=LP[i]+20*i+30)
            k_UP = obtain_peak_linshape(path=path, freq_start=UP[i]-20*i-30, freq_end=UP[i]+20*i+30)
            kLPs.append(k_LP)
            kUPs.append(k_UP)
        # save data
        kLPs, kUPs = np.array(kLPs), np.array(kUPs)
        data = np.zeros((kLPs.size, 4))
        data[:,0] = UP - LP
        data[:,1] = k_bare
        data[:,2] = np.abs(kLPs)
        data[:,3] = np.abs(kUPs)
        np.savetxt(filename, data, fmt='%.4E', header="Rabi splitting, Bare C=O 1/tau, LP 1/tau, UP 1/tau (units cm-1)")
    return UP-LP, k_bare, kLPs, kUPs

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
    #yfit_phx, tau = fit_exponential(t[300:], phx[300:])
    return tau

def obtain_polariton_lifetime_versus_OmegaN(Amp):
    filename = "polariton_tau_2_OmegaN_%s.tmp" %Amp
    try:
        data = np.loadtxt(filename)
        print("Have calculated polariton lifetime versus OmegaN under Amp = %s" %Amp)
        OmegaN, tau_LPs, tau_UPs = data[:,0], data[:,1], data[:,2]
    except:
        print("From stratch fitting polariton lifetime from NE-CavMD under Amp = %s" %Amp)
        omega = 2320
        couplings = ["5e-5", "1e-4", "2e-4", "3e-4", "4e-4", "5e-4"]
        y1 = np.array([2350, 2375, 2428, 2485, 2548, 2616])
        y2 = np.array([2301, 2280, 2242, 2207, 2177, 2150])
        OmegaN = y1 - y2
        paths_LP = ["exciteCO2LP_coupling_%s/Amp_%s" %(c, Amp) for c in couplings]
        paths_LP[2] = "exciteCO2LP/Amp_%s" % Amp
        paths_UP = ["exciteCO2UP_coupling_%s/Amp_%s" %(c, Amp) for c in couplings]
        paths_UP[2] = "exciteCO2UP/Amp_%s" % Amp
        tau_LPs, tau_UPs = [], []
        for i in range(len(couplings)):
            if i == 5 and Amp == "6e-3":
                tau_LP = obtain_polariton_lifetime(path=paths_LP[i], omega=omega, nstart=350)
            else:
                tau_LP = obtain_polariton_lifetime(path=paths_LP[i], omega=omega)
            tau_UP = obtain_polariton_lifetime(path=paths_UP[i], omega=omega)
            tau_LPs.append(tau_LP)
            tau_UPs.append(tau_UP)
        tau_LPs = np.array(tau_LPs)
        tau_UPs = np.array(tau_UPs)
        # save data
        data = np.zeros((tau_LPs.size, 3))
        data[:,0], data[:,1], data[:,2] = OmegaN, tau_LPs, tau_UPs
        np.savetxt(filename, data, fmt='%.4E', header="OmegaN [cm-1], tau_LP [ps], tau_UP [ps]")
    return OmegaN, tau_LPs, tau_UPs

def obtain_polariton_lifetime_versus_detuning(Amp):
    filename = "polariton_tau_2_detuning_%s.tmp" %Amp
    try:
        data = np.loadtxt(filename)
        print("Have calculated polariton lifetime versus OmegaN under Amp = %s" %Amp)
        omegas, tau_LPs, tau_UPs = data[:,0], data[:,1], data[:,2]
    except:
        print("From stratch fitting polariton lifetime from NE-CavMD under Amp = %s" %Amp)
        omegas = [2100, 2150, 2200, 2250, 2300, 2320, 2350, 2400, 2500]
        paths_LP = ["exciteCO2LP_%d/Amp_%s" %(omega, Amp) for omega in omegas]
        paths_LP[5] = "exciteCO2LP/Amp_%s" % Amp
        paths_UP = ["exciteCO2UP_%d/Amp_%s" %(omega, Amp) for omega in omegas]
        paths_UP[5] = "exciteCO2UP/Amp_%s" % Amp
        tau_LPs, tau_UPs = [], []
        for i, omega in enumerate(omegas):
            tau_LP = obtain_polariton_lifetime(path=paths_LP[i], omega=omega)
            tau_UP = obtain_polariton_lifetime(path=paths_UP[i], omega=omega)
            tau_LPs.append(tau_LP)
            tau_UPs.append(tau_UP)
        omegas = np.array([float(x) for x in omegas])
        tau_LPs = np.array(tau_LPs)
        tau_UPs = np.array(tau_UPs)
        # save data
        data = np.zeros((tau_LPs.size, 3))
        data[:,0], data[:,1], data[:,2] = omegas, tau_LPs, tau_UPs
        np.savetxt(filename, data, fmt='%.4E', header="omega_c [cm-1], tau_LP [ps], tau_UP [ps]")
    return omegas, tau_LPs, tau_UPs

'''
def plot_polariton_versus_OmegaN(Amp="6e-4"):
    OmegaN, tau_LPs, tau_UPs = obtain_polariton_lifetime_versus_OmegaN(Amp=Amp)
    xs = [OmegaN]*2
    ys = [1.0/tau_LPs, 1.0/tau_UPs]
    colors = ["r^", "b^"]
    labels = ["LP NE-CavMD", "UP NE-CavMD"]

    # Add plotting for fitting the linear response theory
    OmegaN, k_bare, k_LP, k_UP = collect_polariton_lineshape()
    ys_LR = [k_LP/33.35641, k_UP/33.35641]
    # Adding plotting for 1D model
    theta_half = np.arctan(2.0 * 2320.0/OmegaN) / 2.0
    k_LP = k_bare * np.sin(theta_half)**2.0
    k_UP = k_bare * np.cos(theta_half)**2.0
    ys_1d = [k_LP/33.35641, k_UP/33.35641]

    axes = clp.initialize(2, 1, width=4.3, height=4.3*0.618*2, LaTeX=True, fontsize=13, sharex=True)
    clp.plotone(xs, ys_LR, axes[0], colors=["ro", "bo"], labels=["LP CavMD LR", "UP CavMD LR"])
    clp.plotone(xs, ys_1d, axes[0], colors=["r", "b"], labels=["LP 1D model", "UP 1D model"], alpha=0.6,
            ylabel="decay rate [ps$^{-1}$]")
    clp.plotone(xs, ys, axes[1], labels=labels, colors=colors, xlabel="Rabi splitting [cm$^{-1}$]",
            ylabel="decay rate [ps$^{-1}$]", markersize=10)
    # Also plot a Lorentzian fit for you
    yfit_UP, __ = fit_Lorentzian_centered(UP, 1.0/tau_UPs)
    yfit_LP, __ = fit_biLorentzian(LP, 1.0/tau_LPs)
    clp.plotone(xs, [yfit_LP, yfit_UP], axes[1], colors=["r-.", "b-."], labels=["LP bi-Lorentz fit", "UP Lorentz fit"], alpha=0.6, lw=3)
    from matplotlib.patches import Rectangle
    axes[1].add_patch( Rectangle((OmegaN[2]-25, ys[0][2]-0.25),
                        50, 0.5, color='m', alpha=0.2))
    clp.adjust(savefile="test.pdf")
'''

'''
def plot_polariton_versus_detuning(Amp="6e-4"):
    Omegas, tau_LPs, tau_UPs = obtain_polariton_lifetime_versus_detuning(Amp=Amp)
    xs = [Omegas-2320]*2
    ys = [1.0/tau_LPs, 1.0/tau_UPs]
    colors = ["b-^", "r-o"]
    labels = ["LP strong", "UP strong"]
    ax = clp.initialize(1, 1, width=5., LaTeX=True, fontsize=12)
    clp.plotone(xs, ys, ax, labels=labels, colors=colors, ylog=True, xlabel="cavity mode detuning [cm$^{-1}$]", ylabel="polariton decay rate $1/\\tau$ [ps$^{-1}$]")
    #ax.axvspan(OmegaN[2]-20, OmegaN[2]+20, alpha=0.2, color='c')
    clp.adjust(savefile="test2.pdf")

def plot_photon():
    axes = clp.initialize(1, 2, width=4.3*1.5, height=4.3*0.8, LaTeX=True, fontsize=12, labelthem=True, labelthemPosition=[-0.03, 1.12])
    subpath = "2500/Amp_6e-3"
    paths = ["exciteCO2LP_%s" %subpath, "exciteCO2UP_%s" %subpath]
    #paths = ["exciteCO2LP/%s" %subpath, "exciteCO2UP/%s" %subpath]
    for i in range(2):
        t, phx, phy = prepare_ph_data(path=paths[i], omega=2320 if i is 0 else 2200)
        yfit_phx, tau = fit_exponential(t[300:], phx[300:] + phy[300:])
        xs = [t, t[300:]]
        ys = [phx + phy, yfit_phx]
        colors = ["r", "k-.", "c"]
        labels = ["$x$-polarized ph", "fit for $x$-polarized ph", "$y$-polarized ph"]
        if i == 0:
            clp.plotone(xs, ys, axes[i], colors=colors, labels=labels, xlim=[0.,60], xlabel="time [ps]", ylabel="Photonic potential energy [a.u.]")
            axes[i].text(0.5, 0.5, "$\\tau =$ %.2f ps" %tau, fontsize=12, transform = axes[i].transAxes)
        else:
            clp.plotone(xs, ys, axes[i], colors=colors, showlegend=False, xlim=[0.,60], xlabel="time [ps]")
            axes[i].text(0.5, 0.5, "$\\tau =$ %.2f ps" %tau, fontsize=12, transform = axes[i].transAxes)
        axes[i].axvspan(0.1, 0.6, alpha=0.2, color='y')
    clp.adjust(savefile="photon_dynamics_LP.pdf", tight_layout=True)
    ax = clp.initialize()
    x = np.array([0, 5e-5, 1e-4, 2e-4, 3e-4, 4e-4, 5e-4, 6e-4])
    y1 = np.array([2320, 2350, 2375, 2428, 2485, 2548, 2616, 2690])
    y2 = np.array([2320, 2301, 2280, 2242, 2207, 2177, 2150, 2130])
    xs = [x]
    ys = [y1 - y2]
    clp.plotone(xs, ys, ax, colors=["ro", "k^"])
    clp.adjust()

'''

def plot_lifetime_all():
    #Here I plot 6 sub figures which conlcudes the main findings about polaritonic lifetime
    fig, axes = clp.initialize(1, 2, width=4.3*2, height=4.3*0.618, LaTeX=True, fontsize=12, labelthem=True, labelthemPosition=[-0.03, 1.03],  return_fig_args=True)

    # first plot  #photonic energy as a function of time for UP and LP under weak excitations
    subpath = "Amp_6e-4"
    paths = ["exciteCO2UP/%s" %subpath, "exciteCO2LP/%s" %subpath]
    #paths = ["C13O2CO2_Experiments/exciteC13O2CO2LP_2300/%s" %subpath, "C13O2CO2_Experiments/exciteC13O2CO2LP_2300/%s" %subpath]
    xs, ys, taus = [], [], []
    for i in range(2):
        t, phx, phy = prepare_ph_data(path=paths[i], omega=2320)
        yfit_phx, tau = fit_exponential(t[300:20000], phx[300:20000] + phy[300:20000])
        xs += [t, t[300:20000]]
        ys += [phx + phy, yfit_phx]
        taus.append(tau)
    colors = ["m", "k--", "c", "k-."]
    labels = ["UP NE-CavMD", "UP exp fit", "LP NE-CavMD", "LP exp fit"]
    clp.plotone(xs, ys, axes[0], colors=colors, labels=labels, lw=1.2, xlabel="time [ps]", ylabel="photonic energy [a.u.]", xlim=[0, 5], ylim=[0, 0.0149])
    axes[0].text(0.2, 0.66, "$\\tau_{UP}$ = %.1f ps" %taus[0], fontsize=12, transform=axes[0].transAxes, color=colors[0])
    axes[0].text(0.2, 0.06, "$\\tau_{LP}$ = %.1f ps" %taus[1], fontsize=12, transform=axes[0].transAxes, color=colors[2])
    axes[0].text(0.13, 0.87, E0string_weak, fontsize=12, transform=axes[0].transAxes)
    axes[0].axvspan(0.1, 0.6, alpha=0.1, color='y')

    # second plot  #photonic energy as a function of time for UP and LP under strong excitations

    subpath = "Amp_6e-3"
    paths = ["exciteCO2UP/%s" %subpath, "exciteCO2LP/%s" %subpath]
    #paths = ["C13O2CO2_Experiments/exciteC13O2CO2LP_2300/%s" %subpath, "C13O2CO2_Experiments/exciteC13O2CO2LP_2300/%s" %subpath]
    xs, ys, taus = [], [], []
    for i in range(2):
        t, phx, phy = prepare_ph_data(path=paths[i], omega=2320)
        yfit_phx, tau = fit_exponential(t[300:20000], phx[300:20000] + phy[300:20000])
        xs += [t, t[300:20000]]
        ys += [phx + phy, yfit_phx]
        taus.append(tau)
    colors = ["m", "k--", "c", "k-."]
    labels = ["UP NE-CavMD", "UP exp fit", "LP NE-CavMD", "LP exp fit"]
    clp.plotone(xs, ys, axes[1], colors=colors, labels=labels, lw=1.2, xlabel="time [ps]", ylabel="photonic energy [a.u.]", xlim=[0, 5], ylim=[0, 1.48])
    axes[1].text(0.2, 0.63, "$\\tau_{UP}$ = %.1f ps" %taus[0], fontsize=12, transform=axes[1].transAxes, color=colors[0])
    axes[1].text(0.2, 0.2, "$\\tau_{LP}$ = %.1f ps" %taus[1], fontsize=12, transform=axes[1].transAxes, color=colors[2])
    axes[1].text(0.13, 0.87, E0string_strong, fontsize=12, transform=axes[1].transAxes)
    axes[1].axvspan(0.1, 0.6, alpha=0.1, color='y')
    
    '''
    # Third and forth plot, UP and LP lifetime as a function of Rabi splitting
    OmegaN, tau_LPs_weak, tau_UPs_weak = obtain_polariton_lifetime_versus_OmegaN(Amp="6e-4")
    OmegaN, tau_LPs_strong, tau_UPs_strong = obtain_polariton_lifetime_versus_OmegaN(Amp="6e-3")
    xs = [OmegaN]*6
    # Lorentzian fit
    yfit_UP_weak, __ = fit_Lorentzian_centered(UP, 1.0/tau_UPs_weak)
    yfit_LP_weak, __ = fit_biLorentzian(LP, 1.0/tau_LPs_weak)
    yfit_UP_strong, __ = fit_Lorentzian_centered(UP, 1.0/tau_UPs_strong)
    yfit_LP_strong, __ = fit_biLorentzian(LP, 1.0/tau_LPs_strong)
    ys_LP = [1.0/tau_LPs_weak, 1.0/tau_LPs_strong, yfit_LP_weak,  yfit_LP_strong]
    ys_UP = [1.0/tau_UPs_weak, 1.0/tau_UPs_strong, yfit_UP_weak,  yfit_UP_strong]
    # linear response theory
    OmegaN, k_bare, k_LP, k_UP = collect_polariton_lineshape()
    cminv2psinv = 1.0 / 33.35641
    ys_LP.append(k_LP * cminv2psinv) # cm-1 to ps-1
    ys_UP.append(k_UP * cminv2psinv)
    # 1D model
    theta_half = np.arctan2(2.0 * 2320.0*OmegaN, (OmegaN**2 + 2327**2 - 2320**2)) / 2.0
    k_LP = k_bare * np.sin(theta_half)**2.0
    k_UP = k_bare * np.cos(theta_half)**2.0
    print("cos^2", np.cos(theta_half)**2.0)
    ys_LP.append(k_LP * cminv2psinv)
    ys_UP.append(k_UP * cminv2psinv)

    labels = [E0string_weak,
              E0string_strong, "Lorentzian fit", "Lorentzian fit",
              "CavMD LR", "1D model"]
    colors = ["m*", "mo", "0.5", "0.5", "r^", "g--"]
    clp.plotone(xs[0:3], ys_UP[0:3], axes[1,0], labels=labels, colors=colors, lw=2., markersize=10,
                xlabel="Rabi freq [cm$^{-1}$]", ylabel="UP decay rate [ps$^{-1}$]" , ylog=True, ylim=[0.01, 9.9])
    clp.plotone(xs[3:4], ys_UP[3:4], axes[1,0], labels=labels[3:], colors=colors[3:], lw=2.,
                markersize=7, showlegend=False, alphaspacing=0.3)
    #axes[1,0].text(OmegaN[0], ys_UP[-2][1]*0.3, "CavMD LR", color="r", alpha=0.7)
    #axes[1,0].text(OmegaN[2]+15, ys_UP[-1][-2]*1.1, "1D model", color="g", alpha=0.4)

    colors = ["c*", "co", "0.3", "0.3", "r^", "g--"]
    labels = [E0string_weak,
              E0string_strong, "bi-Lorentzian fit", "bi-Lorentzian fit",
              "CavMD LR", "1D model"]
    clp.plotone(xs[0:3], ys_LP[0:3], axes[1,1], labels=labels, colors=colors, lw=2., markersize=10,
                xlabel="Rabi freq [cm$^{-1}$]", ylabel="LP decay rate [ps$^{-1}$]", ylog=True, ylim=[0.01, 9.9])
    clp.plotone(xs[3:4], ys_LP[3:4], axes[1,1], labels=labels[3:], colors=colors[3:], lw=2.,
                markersize=7, showlegend=False, alphaspacing=0.3)


    ax2 = axes[1,0].twiny()
    ax2.set_xlim(axes[1,0].get_xlim())
    ax2.set_xticks(OmegaN)
    ax2.set_xticklabels(UP, c='m')
    #ax2.set_xlabel("UP freq [cm$^{-1}$]", c="m")
    axes[1,0].text(0.2, 0.9, "UP freq [cm$^{-1}$]", c="m", fontsize=12, transform=axes[1,0].transAxes)

    ax2 = axes[1,1].twiny()
    ax2.set_xlim(axes[1,1].get_xlim())
    ax2.set_xticks(OmegaN)
    ax2.set_xticklabels(LP, c='c')
    #ax2.set_xlabel("LP freq [cm$^{-1}$]", c='c')
    axes[1,1].text(0.4, 0.9, "LP freq [cm$^{-1}$]", c="c", fontsize=12, transform=axes[1,1].transAxes)
    '''
    # New
    '''
    # Next plot for polaritonic decay rate as a function of detuning
    Amps = ["6e-4", "6e-3"]
    Omegas, tau_LPs_weak, tau_UPs_weak = obtain_polariton_lifetime_versus_detuning(Amp="6e-4")
    Omegas, tau_LPs_strong, tau_UPs_strong = obtain_polariton_lifetime_versus_detuning(Amp="6e-3")
    xs = [Omegas-2327]*3

    ys_UP = [1.0/tau_UPs_weak, 1.0/tau_UPs_strong]
    ys_LP = [1.0/tau_LPs_weak, 1.0/tau_LPs_strong]

    # 1D model
    print("cavity frequency", Omegas)
    theta_half = np.arctan2(2.0 * Omegas*OmegaN[2], (OmegaN[2]**2 + 2327**2 - Omegas**2)) / 2.0
    k_LP = k_bare * np.sin(theta_half)**2.0
    k_UP = k_bare * np.cos(theta_half)**2.0
    print("sin^2", np.sin(theta_half)**2.0)
    ys_LP.append(k_LP * cminv2psinv)
    ys_UP.append(k_UP * cminv2psinv)

    colors = ["m*", "mo", "g--"]
    labels = [E0string_weak, E0string_strong, "1D model"]
    clp.plotone(xs, ys_UP, axes[2,0], labels=labels, colors=colors, markersize=10, ylim=[7e-3, 90],
            ylog=True, xlabel="cavity mode detuning [cm$^{-1}$]", ylabel="UP decay rate [ps$^{-1}$]",
            alphaspacing=0.2)

    colors = ["c*", "co", "g--"]
    labels = [E0string_weak, E0string_strong, "1D model"]
    clp.plotone(xs, ys_LP, axes[2,1], labels=labels, colors=colors, markersize=10, ylim=[7e-3, 90],
            ylog=True, xlabel="cavity mode detuning [cm$^{-1}$]", ylabel="LP decay rate [ps$^{-1}$]",
            alphaspacing=0.2)
    '''
    
    clp.adjust(savefile="lifetime_all_half.pdf", tight_layout=True)

if __name__ == "__main__":
    plot_lifetime_all()
