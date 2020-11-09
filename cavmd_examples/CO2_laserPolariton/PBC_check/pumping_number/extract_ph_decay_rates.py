import numpy as np
import glob
import sys
import columnplots as clp
from scipy.optimize import curve_fit

def func(x, k, a, b):
    return a * np.exp(-x/k) + b

def fit_exponential(x, y):
    popt, pocv =  curve_fit(func, x, y)
    yfit = func(x, *popt)
    print("popt is", popt)
    return yfit, popt[0]

def prepare_ph_data(path="pumping_Rabi/E0_2e-4_Exc_UP_Amp_6e-3", omega=2320, pattern="simu_*.out", dtfs=2):
    filenames = glob.glob("%s/%s" %(path, pattern))
    N = len(filenames)
    print("reading %s with %d files" %(path, N))
    size = np.size(np.loadtxt(filenames[0])[:,0])
    t, phx, phy = np.zeros(size), np.zeros(size), np.zeros(size)
    ph2au = 0.5 * 1.8897259885789**2 * omega**2 * 0.0000045563352812122295**2
    vph2au = 0.5 * 1.8897259885789**2 / 41.341374575751**2
    for filename in filenames:
        data = np.loadtxt(filename)
        t_loc, qx_loc, qy_loc = data[:,1], data[:,-6], data[:,-2]
        vx_loc = np.gradient(qx_loc, dtfs, axis=-1, edge_order=2)
        vy_loc = np.gradient(qy_loc, dtfs, axis=-1, edge_order=2)
        t += t_loc
        phx += qx_loc**2 * ph2au + vx_loc**2 * vph2au
        phy += qy_loc**2 * ph2au + vy_loc**2 * vph2au
    return t/N, (phx + phy)/N

def extract_ph_decay_rates(paths, omega_lst, PlotIt=False):
    ks = []
    for omega, path in zip(omega_lst, paths):
        # check if the decay rate has been calculated
        try:
            data_tau = np.loadtxt(path+"/ph_lifetime.txt")
            tau = data_tau
            print("tau = %.2f ps reading from saved data" %tau)
        except:
            t, ph = prepare_ph_data(path=path, omega=omega)
            yfit_ph, tau = fit_exponential(t[300:], ph[300:])
            print("tau = %.2f ps reading from calculating" %tau)
            data_tau = np.array([tau])
            np.savetxt(path+"/ph_lifetime.txt", data_tau)
        if PlotIt:
            plot_single_ph_traj(path=path, omega=omega)
        ks.append(1.0/tau)
    ks = np.array(ks)
    return ks

def plot_single_ph_traj(path="pumping_Rabi/E0_2e-4_Exc_UP_Amp_6e-3", omega=2320):
    t, ph = prepare_ph_data(path=path, omega=omega)
    ax = clp.initialize(1, 1)
    clp.plotone([t], [ph], ax, showlegend=False)
    clp.adjust()

if __name__ == '__main__':
    paths = sys.argv[1:]
    extract_ph_decay_rates(paths=paths, omega_lst=[2320]*len(paths), PlotIt=False)
