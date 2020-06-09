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
    paths = ["E0_0e-4", "E0_2e-4", "E0_4e-4", "E0_6e-4", "E0_8e-4"]
    #paths = ["E0_0e-4", "E0_1e-4", "E0_3e-4", "E0_5e-4", "E0_7e-4"]
    colors = ["k", "g", "c", "r", "b"]
    labels = ["outside cavity", "in cavity $\widetilde{\\varepsilon}=2\\times 10^{-4}$",
                "in cavity $\widetilde{\\varepsilon}=4\\times 10^{-4}$",
                "in cavity $\widetilde{\\varepsilon}=6\\times 10^{-4}$",
                "in cavity $\widetilde{\\varepsilon}=8\\times 10^{-4}$"]
    axes = clp.initialize(5, 1, width=4.3, height=6*1.25, LaTeX=True, fontsize=12, sharex=True,
            labelthem=True, labelthemPosition=[0.98, 0.92],
            commonX=[0.5, 0.05, "frequency [cm$^{-1}$]"],
            commonY=[0.03, 0.5, r"$n(\omega)\alpha(\omega)$ [arb. units]"])
    LPs = [3424, 3304, 3229, 3150]
    UPs = [3844, 4039, 4339, 4661]
    for i, path in enumerate(paths):
        data = obtain_avg_data(path, pattern="simu_*.dac.txt")
        x, y = data[:,5], (data[:,6] + data[:,7])/2e28
        if i is 0 or i is 1:
            ylim = [-0.035, 1.62]
        else:
            ylim = None
        clp.plotone([x], [y], axes[i], colors=[colors[i]], labels=[labels[i]], lw=1.5, xlim=[0, 5300], ylim=ylim)
        axes[i].legend(loc="upper left", bbox_to_anchor=(0.01,0.99))
        if i > 0:
            axes[i].text(LPs[i-1]-10, 1.0, "LP")
            axes[i].text(UPs[i-1]+60, i, "UP")
    clp.adjust(savefile="IR.pdf")

if __name__ == "__main__":
    plot_IR()
