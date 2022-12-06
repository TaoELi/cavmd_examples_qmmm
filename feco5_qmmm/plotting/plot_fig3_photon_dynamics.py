import numpy as np
import columnplots as clp
import MDAnalysis as mda
import pandas as pd
from scipy import signal
import os, sys, time
import math
from itertools import islice
from scipy import fftpack
import glob
import json
import matplotlib.cm as cm
import matplotlib.transforms as mtransforms
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit

au2cminv = 219474.63
e_KT = 208.5
photon_energy_cminv = 2018.0

def obtain_avg_data(path, pattern="simu_*.vac.txt"):
    filenames = glob.glob("%s/%s" %(path, pattern))
    print("Averaging over %d trajectories under %s" %(len(filenames), path))
    data = np.loadtxt(filenames[0])
    for filename in filenames[1:]:
        data += np.loadtxt(filename)
    data /= float(len(filenames))
    return data

def func(x, k, a, b):
    return a * np.exp(-x/k) + b

def func_biexp(x, tau1, tau2, a, b):
    return -a * np.exp(x/tau1) + b * np.exp(-x/tau2)

def fit_exponential(x, y):
    popt, pocv =  curve_fit(func, x, y)
    yfit = func(x, *popt)
    print("popt is", popt)
    return yfit

def fit_biexponential(x, y):
    popt, pocv =  curve_fit(func_biexp, x, y)
    yfit = func_biexp(x, *popt)
    print("popt is", popt)
    return yfit

def construct_velo_qmmm(data, idx, N=16):
    data_slice = data[idx,:]
    velo = np.zeros((N*11+2, 3))
    # get data of molecules
    data_molec = data_slice[18:]
    #print("molecular data have size of", data_molec.size)
    data_molec = np.reshape(data_molec, (N*11, 3))
    data_ph = data_slice[12:18]
    data_ph = np.reshape(data_ph, (2, 3))
    velo[0:N*11,:] = data_molec
    velo[N*11:,:] = data_ph
    return velo

def get_ph_energy(velo):
    v_q1 = velo[-2,:]
    v_q2 = velo[-1,:]
    e = 0.5 * np.sum(v_q1**2 + v_q2**2)
    return e * au2cminv

def moving_average(arr, window_size=20):
    # Convert array of integers to pandas series
    numbers_series = pd.Series(arr)
    # Get the window of series
    # of observations of specified window size
    windows = numbers_series.rolling(window_size, center=True)
    # Create a series of moving
    # averages of each window
    moving_averages = windows.mean()
    # Convert pandas series back to list
    moving_averages_list = np.array(moving_averages.tolist())
    return moving_averages_list

def get_qmmm_ph_data(filename, N=16):
    data = np.loadtxt(filename)
    nframes = np.shape(data)[0]
    print("There are %d qmmm frames" %nframes)
    e_ph = []
    for idx in range(nframes):
        current_veloc = construct_velo_qmmm(data, idx, N=N)
        e_ph.append(get_ph_energy(current_veloc))
    e_ph = np.array(e_ph) - e_KT*3
    t =  np.array([2.0*i for i in range(e_ph.size)])
    return t*1e-3, e_ph

def plot_photon_dynamics_qmmm(window_size=20):
    fig, axes = clp.initialize(1, 2, width=4.3*2, height=4.3*0.618, LaTeX=True, fontsize=12, labelthem=True, labelthemPosition=[-0.03, 1.03],  return_fig_args=True)

    # first plot  #photonic energy as a function of time for UP and LP under weak excitations
    subpath = "../qmmm/incav_eq_weak_bath"
    paths = ["%s/E0_1e-4_excite_UP_weaker/" %subpath, "%s/E0_1e-4_excite_LP_weaker/" %subpath]
    xs, ys, xfits, yfits = [], [], [], []
    for i in range(2):
        t, e_ph = get_qmmm_ph_data(filename=paths[i]+"simu_1.out")
        e_ph = moving_average(e_ph, window_size=window_size)
        # truncate NaN
        t = t[~np.isnan(e_ph)]
        e_ph = e_ph[~np.isnan(e_ph)]
        print(np.shape(t))
        yfit = fit_exponential(t[500:], e_ph[500:] / photon_energy_cminv)
        yfits.append(yfit)
        xfits.append(t[500:])
        xs += [t]
        ys += [e_ph / photon_energy_cminv]
    colors = ["m",  "c",]
    labels = ["UP", "LP"]
    clp.plotone(xs, ys, axes[0], colors=colors, labels=labels, lw=1.2, xlabel="time [ps]", ylabel="photon KE [$\hbar\omega_c$]", xlim=[0, 2.5], ylim=[0, 17])
    #clp.plotone(xfits, yfits, axes[0], colors=["k--", "k--"], lw=1, showlegend=False)
    E0string_weak = "$n$-dodecane 300K\npluse amp. $E_0=1.7\\times 10^{-3}$ a.u."
    axes[0].text(0.05, 0.8, E0string_weak, fontsize=12, transform=axes[0].transAxes)
    axes[0].axvspan(0., 1., alpha=0.1, color='y')

    # second plot  #photonic energy as a function of time for UP and LP under weak excitations
    subpath = "../qmmm/incav_eq_weak_bath"
    paths = ["%s/E0_1e-4_excite_UP_stronger/" %subpath, "%s/E0_1e-4_excite_LP_stronger/" %subpath]
    xs, ys, xfits, yfits = [], [], [], []
    for i in range(2):
        t, e_ph = get_qmmm_ph_data(filename=paths[i]+"simu_1.out")
        e_ph = moving_average(e_ph, window_size=window_size)
        # truncate NaN
        t = t[~np.isnan(e_ph)]
        e_ph = e_ph[~np.isnan(e_ph)]
        print(np.shape(t))
        yfit = fit_exponential(t[500:], e_ph[500:] / photon_energy_cminv)
        yfits.append(yfit)
        xfits.append(t[500:])
        xs += [t]
        ys += [e_ph / photon_energy_cminv]
    colors = ["m",  "c",]
    labels = ["UP", "LP"]
    clp.plotone(xs, ys, axes[1], colors=colors, labels=labels, lw=1.2, xlabel="time [ps]", xlim=[0, 2.5], showlegend=False, ylim=[0,68])
    #clp.plotone(xfits, yfits, axes[1], colors=["k--", "k--"], lw=1, showlegend=False)
    E0string_strong = "$n$-dodecane 300K\npluse amp. $E_0=3.4\\times 10^{-3}$ a.u."
    axes[1].text(0.05, 0.8, E0string_strong, fontsize=12, transform=axes[1].transAxes)
    axes[1].axvspan(0., 1., alpha=0.1, color='y')

    clp.adjust(savefile="photon_dynamics.pdf", tight_layout=True)

if __name__ == "__main__":
    plot_photon_dynamics_qmmm()
