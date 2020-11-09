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


def plot_diffusion():
    paths = ["E0_0e-4", "E0_4e-4"]
    colors = ["k", "c--"]
    labels = ["outside cavity", "in cavity $\widetilde{\\varepsilon}=4\\times 10^{-4}$"]
    axes = clp.initialize(2, 1, width=4.3, height=5, LaTeX=True, fontsize=12,
            labelthem=True, labelthemPosition=[0.98, 0.92])
    ts, corrs, freqs, sps = [], [], [], []
    for i, path in enumerate(paths):
        data = obtain_avg_data(path, pattern="simu_*.vac.txt")
        t, corr, freq, sp = data[:,0], data[:, 4], data[:,5], data[:,-1]
        ts.append(t*1e-3)
        corrs.append(corr*1e6)
        print("diffusion constatn is %.3E" %(np.sum(corr) * (t[2]-t[1]) * 1e3 / 3.0))
        freqs.append(freq)
        sps.append(sp*2e3)

    clp.plotone(ts, corrs, axes[0], colors=colors, labels=labels, lw=1.5,
        xlabel="time [ps]", ylabel=r"$C_{vv}(t)$ [$\AA^2$/ps$^2$]", xlim=[-0.02, 1])
    clp.plotone(freqs, sps, axes[1], colors=colors, labels=labels, lw=1.5, xlim=[-100, 5000],
        xlabel="frequency [cm$^{-1}$]", ylabel=r"$C_{vv}(\omega)$ [arb. units]", showlegend=False)
    axes[0].legend(loc="upper left", bbox_to_anchor=(0.1,0.99))
    clp.adjust(savefile="Diffusion.pdf", tight_layout=True)

if __name__ == "__main__":
    plot_diffusion()
