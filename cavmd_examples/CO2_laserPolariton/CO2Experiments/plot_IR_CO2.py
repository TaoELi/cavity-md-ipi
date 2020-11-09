import numpy as np
import columnplots as clp
from scipy import signal
import os, sys, time
import math
from itertools import islice
from scipy import fftpack
import glob
import json
import matplotlib.cm as cm
import matplotlib.pyplot as plt


def obtain_avg_data(path, pattern="simu_*.vac.txt"):
    filenames = glob.glob("%s/%s" %(path, pattern))
    print("Averaging over %d trajectories under %s" %(len(filenames), path))
    data = np.loadtxt(filenames[0])
    for filename in filenames[1:]:
        data += np.loadtxt(filename)
    data /= float(len(filenames))
    return data

def obtain_IR(path):
    data = obtain_avg_data(path=path, pattern="simu_*.dac.txt")
    freq = data[:,5]
    sp = (data[:,6] + data[:,7]) / 2e28
    return freq, sp

def plot_IR_all():
    '''
    A function to plot Setup, Rabi splitting, and Avoid crossing at the same time
    '''
    fig, axes = clp.initialize(1, 3, width=4.3*3, height=4.3*0.7,  LaTeX=True, fontsize=12, return_fig_args=True, labelthem=True, labelthemPosition=[-0.03, 1.06])

    # First figure, reserve for setup demonstration
    axes[0].axis('off')
    im = plt.imread('setup.png')
    nx, ny, d = np.shape(im)
    dnx = int(nx * 0.02)//2
    dny = int(nx * 0.52)//2
    axes[0].imshow(im[dnx:-dnx,dny:-dny,:])

    # Second figure, plot Rabi splitting as a function of effective coupling strength
    E0s = ["0e-4", "5e-5", "1e-4", "2e-4", "3e-4", "4e-4", "5e-4"]
    paths = ["CO2only_changeE0/E0_%s" %E0 for E0 in E0s]
    xs, ys = [], []
    for i, path in enumerate(paths):
        freq, sp = obtain_IR(path=path)
        xs.append(freq)
        ys.append(sp +  i*2)
    colors = ["k"] * len(path)
    clp.plotone(xs, ys, axes[1], colors=colors, showlegend=False, lw=1, xlim=[1900, 2650], ylim=[-0.2, 24.],
                xlabel="freq [cm$^{-1}$]", ylabel="$n(\omega)\\alpha(\omega)$ [arb. units]")
    colormap = plt.cm.hot
    colors = [colormap(i) for i in np.linspace(0, 0.6,len(axes[1].lines))]
    for i,j in enumerate(axes[1].lines):
        j.set_color(colors[i])

    axes[1].text(1910, 0.3, "outside cavity", fontsize=8)
    axes[1].text(1910, 2.4, "$\widetilde{\\varepsilon} = 5\\times 10^{-5}$ a.u.",fontsize=8)
    for i in range(len(paths[2:])):
        axes[1].text(1910, i*2 + 4.4, "$\widetilde{\\varepsilon} = %d\\times 10^{-4}$ a.u." %(i+1), fontsize=8)

    # identifying Rabi splitting
    df = freq[2] - freq[1]
    nstart, nmiddle, nend = int(2000 // df), int(2327 // df), int(2650 // df)
    print(freq[nstart], freq[nmiddle], freq[nend])
    subfreq1 = freq[nstart:nmiddle]
    subfreq2 = freq[nmiddle:nend]
    OmegaN = [0.0]
    for sp in ys[1:]:
        omega_LP = subfreq1[np.argmax(sp[nstart:nmiddle])]
        omega_UP = subfreq2[np.argmax(sp[nmiddle:nend])]
        print("LP = %.1f, UP = %.1f, Rabi = %.1f" %(omega_LP, omega_UP, omega_UP - omega_LP))
        OmegaN.append(omega_UP - omega_LP)
    E0s = np.array([float(E0) for E0 in E0s])
    OmegaN = np.array(OmegaN)
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    inset_ax2 = inset_axes(axes[1], width="66%", height=0.7)
    clp.plotone([E0s*1e4], [OmegaN], inset_ax2, colors=["r-o"], lw=1.5, showlegend=False)
    inset_ax2.set_xlabel('$\widetilde{\\varepsilon}$ [$\\times 10^{-4}$ a.u.]', fontsize=8)
    inset_ax2.set_ylabel('$\Omega_N$ [cm$^{-1}$]', fontsize=8)
    inset_ax2.tick_params(axis="x", labelsize=8)
    inset_ax2.tick_params(axis="y", labelsize=8)


    # Third figure, plot previous results...
    N = 23
    freqs = [2000 + 25*i for i in range(N)]
    paths = ["CO2only_changeFreq/Freq_%d" %freq for freq in freqs]
    data = obtain_avg_data(path=paths[0], pattern="simu_*.dac.txt")
    nsize = np.shape(data)[0]
    data_all = np.zeros((nsize, N))
    sps = []
    for i, path in enumerate(paths):
        data = obtain_avg_data(path=path, pattern="simu_*.dac.txt")
        data_all[:,i] = (data[:,6] + data[:,7]) / 2e28
        sps.append(data_all[:,i])
    LPs, UPs = [], []
    for freq, sp in zip(freqs, sps):
        omega_LP = subfreq1[np.argmax(sp[nstart:nmiddle])]
        omega_UP = subfreq2[np.argmax(sp[nmiddle:nend])]
        print("Freq = %d, LP = %.1f, UP = %.1f, Rabi = %.1f" %(freq, omega_LP, omega_UP, omega_UP - omega_LP))
        LPs.append(omega_LP)
        UPs.append(omega_UP)
    LPs = np.array(LPs)
    UPs = np.array(UPs)

    x = data[:,5]
    dx = x[2] - x[1]
    nstart, nend = int(2000 / dx), int(2600 / dx)
    x = x[nstart:nend]
    data_all = data_all[nstart:nend,:]

    data_all = np.abs(data_all[::-1, :])
    extent = [freqs[0], freqs[-1], x[0] , x[-1]]

    from matplotlib.colors import LogNorm
    vmax = np.max(np.max(data_all))
    vmin = vmax * 0.02
    pos = axes[2].imshow(data_all, aspect='auto', extent=extent,  cmap=cm.hot,
            interpolation=None,
            norm=LogNorm(vmin=vmin, vmax=vmax))

    xs = [np.array(freqs)]*2
    ys = [np.ones(len(freqs)) * 2327, np.array(freqs)]
    clp.plotone(xs, ys, axes[2], showlegend=False, colors=["w--", "g--"], lw=1.2,
            xlabel="cavity mode freq [cm$^{-1}$]",
            ylabel="IR freq [cm$^{-1}$]")
    #clp.plotone(xs, [LPs, UPs], axes[2], colors=["bo", "bo"], showlegend=False)
    axes[2].text(2030, 2195, "cavity photon", color='g', fontsize=12)
    axes[2].text(2030, 2265, "C=O asym. stretch", color='w', fontsize=12)
    axes[2].tick_params(color='c', labelsize='medium', width=2)

    fig.subplots_adjust(right=0.95)
    cbar_ax = fig.add_axes([0.96, 0.15, 0.02, 0.70])

    cbar = fig.colorbar(pos, cax=cbar_ax)
    cbar.set_label('IR intensity [arb. units]')

    clp.adjust(savefile="IR_CO2.pdf")


if __name__ == "__main__":
    plot_IR_all()
