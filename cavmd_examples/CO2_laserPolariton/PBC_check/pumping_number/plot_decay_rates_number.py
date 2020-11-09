import numpy as np
import glob
import sys
import columnplots as clp
from scipy.optimize import curve_fit
from spectral_overlap import spectral_overlap
from extract_ph_decay_rates import extract_ph_decay_rates

'''
Under weak excitations
Fig. a plots IR spectrum under different effective coupling strengths
Fig. b plots (i) spectral overlap integral and (ii) polariton decay rates versus Rabi splitting
'''

def func(x, k, a):
    return k*x + a

def fit_linear(x, y):
    popt, pocv =  curve_fit(func, x, y)
    xfit = np.concatenate((x, np.array([0.0])))
    yfit = func(xfit, *popt)
    return xfit, yfit

def plot_decay_rates_number():
    N_lst = np.array([216., 280., 343., 432., 625., 800., 1200.])

    def plot_sub(Amp="6e-3", Polariton="LP", ax=None, ylim=None, xlim=None):
        paths = ["N_%d_Exc_%s_Amp_%s/" %(N, Polariton, Amp) for N in N_lst]
        ks = extract_ph_decay_rates(paths=paths, omega_lst=[2320]*len(paths))
        xfit, ks_fit = fit_linear(1.0/N_lst, ks)
        if Amp == "6e-3":
            label="$F = 632$ mJ/cm$^{2}$ \n excite %s" %Polariton
        else:
            label="$F = 6.32$ mJ/cm$^{2}$ \n excite %s" %Polariton
        color = "c" if Polariton == "LP" else "m"
        clp.plotone([1.0/N_lst], [ks], ax, colors=[color+"o"], showlegend=False, ylim=ylim, xlim=xlim)
        clp.plotone([xfit], [ks_fit], ax, colors=["k--"], showlegend=False, lw=1)
        ax.text(0.3, 0.65, label, transform=ax.transAxes, color=color)

    # Define the plot framework
    axes = clp.initialize(2, 2, width=5.3,  fontsize=12, LaTeX=True,
            labelthem=True, labelthemPosition=[0.15, 0.95], sharex=True,
            commonX=[0.5, -0.05, "$1/N_{sub}$"],
            commonY=[-0.05, 0.5, "polariton decay rate [ps$^{-1}$]"])

    plot_sub(Amp="6e-4", Polariton="LP", ax=axes[0,0], xlim=[0, 5e-3], ylim=[0, 2.])
    plot_sub(Amp="6e-4", Polariton="UP", ax=axes[0,1], xlim=[0, 5e-3], ylim=[0, 2.])
    plot_sub(Amp="6e-3", Polariton="LP", ax=axes[1,0], xlim=[0, 5e-3], ylim=[0, 12.])
    plot_sub(Amp="6e-3", Polariton="UP", ax=axes[1,1], xlim=[0, 5e-3], ylim=[0, 2.])

    clp.adjust(tight_layout=True, savefile="decay_rates_number.pdf")


if __name__ == '__main__':
    plot_decay_rates_number()
