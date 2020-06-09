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


def plot_orientation():
    paths = ["E0_0e-4", "E0_4e-4", "E0_6e-4", "E0_8e-4"]
    colors = ["k--", "c", "r--", "b-."]
    labels = ["outside cavity", "in cavity $\widetilde{\\varepsilon}=4\\times 10^{-4}$",
                "in cavity $\widetilde{\\varepsilon}=6\\times 10^{-4}$",
                "in cavity $\widetilde{\\varepsilon}=8\\times 10^{-4}$"]
    fig, axes = clp.initialize(2, 1, width=4.3, height=5, LaTeX=True, fontsize=12,
            labelthem=True, labelthemPosition=[0.1, 0.95], return_fig_args=True)
    ts, corrs, freqs, sps = [], [], [], []
    for i, path in enumerate(paths):
        data = obtain_avg_data(path, pattern="simu_*.oac1.txt")
        t, corr, freq, sp = data[:,0], data[:, 4], data[:,5], data[:,-1]
        #t, corr, freq, sp = data[:,0], data[:, 1], data[:,2], data[:,-1]
        ts.append(t*1e-3)
        corrs.append(corr)
        print("orientation constat is %.3E" %(np.sum(corr) * (t[2]-t[1]) * 1e3))
        freqs.append(freq)
        sps.append(sp / 1e26)

    clp.plotone(ts, corrs, axes[0], colors=colors, labels=labels, lw=1.5,
        xlabel="time [ps]", ylabel=r"$C_{1}^z(t)$ ", xlim=[0, 10], showlegend=False)
    clp.plotone(freqs, sps, axes[1], colors=colors, labels=labels, lw=1.5, xlim=[0, 5000],
        xlabel="frequency [cm$^{-1}$]", ylabel=r"$I_{1}^z(\omega)$ [arb. units]", showlegend=False)
    from mpl_toolkits.axes_grid.inset_locator import inset_axes

    inset_ax1 = inset_axes(axes[0], width="30%", height=1.0, loc=1)
    clp.plotone(ts, corrs, inset_ax1, colors=colors, labels=labels, lw=1.5, xlim=[0, 0.12], ylim=[0.9, 1.01], showlegend=False)

    inset_ax2 = inset_axes(axes[1], width="30%", height=1.0, loc="upper center")
    clp.plotone(freqs, sps, inset_ax2, colors=colors, labels=labels, lw=1.5, xlim=[3890, 4780], ylim=[0., 0.155], showlegend=False)\

    lgd = axes[1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.24),
          fancybox=False, shadow=True, ncol=2)

    clp.adjust(savefile="Orientation.pdf", tight_layout=True, includelegend=lgd)

if __name__ == "__main__":
    plot_orientation()
