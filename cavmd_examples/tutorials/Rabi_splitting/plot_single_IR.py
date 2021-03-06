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
    labels = ["outside cavity", "inside cavity"]
    colors = ["k", "r"]
    ax = clp.initialize(1, 1, width=4.3, LaTeX=False, fontsize=12)
    xs, ys = [], []
    patterns = ["water_*.dac.txt", "simu_*.dac.txt"]
    for pattern in patterns:
        data = obtain_avg_data(path="./", pattern=pattern)
        x, y = data[:,5], (data[:,6] + data[:,7])/2e28
        xs.append(x)
        ys.append(y)
    clp.plotone(xs, ys, ax, colors=colors, labels=labels, lw=1.5, xlim=[0, 5300], xlabel="freq [cm-1]", ylabel="Intensity [arb. units]")
    clp.adjust(savefile="IR_single.png")

if __name__ == "__main__":
    plot_IR()
