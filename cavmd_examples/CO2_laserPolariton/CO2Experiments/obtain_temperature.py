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
    data = np.loadtxt(filenames[1])
    data = np.zeros(np.shape(data))
    N = 0
    for filename in filenames:
        try:
            data += np.loadtxt(filename)
            N += 1
        except:
            print("unfinished simulation for %s" % filename)    
    data /= N
    return data


def plot_temperature(path):
    data = obtain_avg_data(path=path, pattern="simu_*.out")
    t, temp = data[:, 1], data[:, 3]
    T = np.mean(temp[np.size(temp)//2:])
    print("Temperature is %.2f" %T)
    ax = clp.initialize()
    clp.plotone([t], [temp], ax)
    clp.adjust()

if __name__ == "__main__":
    plot_temperature(sys.argv[-1])
