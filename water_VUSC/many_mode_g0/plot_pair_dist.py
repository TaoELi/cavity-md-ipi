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


def plot_pair_distribution():
    paths = ["E0_0e-4", "E0_4e-4", "E0_6e-4"]
    colors = ["k", "c--", "r-."]
    labels = ["outside cavity", "in cavity $\widetilde{\\varepsilon}=4\\times 10^{-4}$", "in cavity $\widetilde{\\varepsilon}=8\\times 10^{-4}$"]
    xs, ys = [], []
    for path in paths:
        data = obtain_avg_data(path, pattern="simu_*.x_*.xyz.pair_dist.txt")
        x, y = data[:,0], data[:,1]
        xs.append(x)
        ys.append(y)

    ax = clp.initialize(1, 1, width=4.3, LaTeX=True, fontsize=12)
    clp.plotone(xs, ys, ax, colors=colors, labels=labels, lw=1.5, xlim=[0, 9],
                xlabel=r"$r_{OO}$ ($\AA$)",
                ylabel="$g(r)$")
    clp.adjust(tight_layout=True, savefile="pair_dist.pdf")

if __name__ == "__main__":
    plot_pair_distribution()
