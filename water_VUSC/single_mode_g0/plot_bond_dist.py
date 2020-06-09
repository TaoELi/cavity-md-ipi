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


def plot_bond_distribution():
    paths = ["E0_0e-4", "E0_4e-4"]
    colors = ["k", "c--"]
    labels = ["outside cavity", "in cavity $\widetilde{\\varepsilon}=4\\times 10^{-4}$"]
    xs, ys = [], []
    for path in paths:
        data = obtain_avg_data(path, pattern="simu_*.bond_length_dist.txt")
        x, y = data[:,0], data[:,1]
        y /= np.sum(y) * (x[2] - x[1])
        xs.append(x)
        ys.append(y)

    ax = clp.initialize(1, 1, width=4.3, LaTeX=True, fontsize=12)
    clp.plotone(xs, ys, ax, colors=colors, labels=labels, lw=1.5, xlim=[0.8, 1.4],
                xlabel=r"O-H bond length ($\AA$)",
                ylabel="Normalized distribution")
    clp.adjust(tight_layout=True, savefile="bond_dist.pdf")

if __name__ == "__main__":
    plot_bond_distribution()
