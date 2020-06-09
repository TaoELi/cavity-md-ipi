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


def plot_IR():
    E0s = [0, 2, 4, 6]
    paths = ["E0_%de-4" %x for x in E0s]
    colors = ["k", "k", "k", "k"]
    labels = ["outside cavity"] + [r"$\widetilde{\varepsilon} = %d \times 10^{-4}$" %x for x in E0s[1:]]
    axes = clp.initialize(4, 1, width=4.3, height=4*1.25, LaTeX=True, fontsize=12, sharex=True,
            labelthem=True, labelthemPosition=[0.10, 0.92],
            commonX=[0.5, 0.05, "frequency [cm$^{-1}$]"],
            commonY=[0.03, 0.5, r"$n(\omega)\alpha(\omega)$ [arb. units]"])
    for i, path in enumerate(paths):
        data = obtain_avg_data(path, pattern="simu_*.dac.txt")
        x, y = data[:,5], (data[:,6] + data[:,7])/2e28
        ylim = [-0.15, 3]
        clp.plotone([x], [y], axes[i], colors=[colors[i]], labels=[labels[i]], lw=1., xlim=[0, 8000], ylim=ylim)
        for j in range(4):
            axes[i].axvline(x=3550.0/4.0*2.0*(j+1), alpha=0.5)
        axes[i].legend(loc="upper right")
    clp.adjust(savefile="IR_multimode_g0.pdf")

if __name__ == "__main__":
    plot_IR()
