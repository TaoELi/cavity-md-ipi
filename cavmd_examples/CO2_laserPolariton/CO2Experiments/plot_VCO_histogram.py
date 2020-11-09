import numpy as np
import columnplots as clp
from scipy import signal
import os, sys, time
import math
import glob
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
from matplotlib import animation


def obtain_tot_data(path, pattern="simu_*.vac.txt"):
    filenames = glob.glob("%s/%s" %(path, pattern))
    print("Averaging over %d trajectories under %s" %(len(filenames), path))
    tot_data = []
    for filename in filenames:
        data = np.loadtxt(filename)
        tot_data.append(data)
    tot_data = np.array(tot_data)
    print(np.shape(tot_data))
    return tot_data

def plot_hist():
    paths = ["exciteCO2LP/Amp_6e-3", "exciteCO2LP_2200/Amp_6e-3", "exciteCO2LP_nocavity/Amp_6e-3"]

    data_tot = []
    for i, path in enumerate(paths):
        data = obtain_tot_data(path=path, pattern="simu_*.xc.xyz.VCO_td.txt") / 0.010603
        data_tot.append(data)

    labels = ["excite LP in 2320 cm$^{-1}$ cavity", "excite LP in 2200 cm$^{-1}$ cavity", "excite 2241 cm$^{-1}$ outside cavity"]
    for h in [20]:
        axes = clp.initialize(col=3, row=1, width=4.3, height=4.3*0.618*2, fontsize=12,
            sharex=True, sharey=True, LaTeX=True, labelthem=True, labelthemPosition=[-0.1, 1.1])
        for i, path in enumerate(paths):
            data = np.reshape(data_tot[i][:,h,1:], -1)
            print(np.shape(data))
            axes[i].hist(data, range=(0, 6), bins=30, density=True, log=True, color = "orchid", ec="black")
            axes[i].set_ylim(1e-3, 10)
            axes[i].set_ylabel("Density distr.")
            axes[i].text(0.3, 0.8, labels[i], fontsize=12, transform=axes[i].transAxes)
        axes[-1].set_xlabel("C=O bond potential energy in CO$_2$ [$\hbar\omega_0$]")
        clp.adjust(tight_layout=True, savefile="VCO_dist_LP.pdf")


if __name__ == "__main__":
    plot_hist()
