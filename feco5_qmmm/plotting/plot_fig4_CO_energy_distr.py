'''
This script collects the CO bond energy distribution at 1.5 ps at three conditions
1) No pulse excitation
2) LP excitation
3) UP excitation
Compare the CO bond energy distribution under these three conditions

Moreover, we also plot CO bond at axial or equatorial positions distribution

Showing energy transfer between localized CO bonds at different locations
'''

import numpy as np
import MDAnalysis as mda
import columnplots as clp
import pandas as pd
from scipy import signal
import os, sys, time
import math
import glob
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
from matplotlib import animation
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

def find_axial_Cs(xyz_feco5, which_mole=0):
    xyz_Fe = xyz_feco5[0,:]
    angle_dict = {}
    for idx_C1 in range(1, 5):
        for idx_C2 in range(idx_C1+1, 6):
            xyz_C1 = xyz_feco5[idx_C1,:]
            xyz_C2 = xyz_feco5[idx_C2,:]
            # calculate the angle C1-Fe-C2
            v1 = xyz_C1 - xyz_Fe
            v2 = xyz_C2 - xyz_Fe
            v1 /= np.linalg.norm(v1)
            v2 /= np.linalg.norm(v2)
            dot_product = np.dot(v1, v2)
            angle1 = np.arccos(dot_product) * 180 / np.pi
            #print("Angle between C%d-Fe-C%d is %.2f" %(idx_C1, idx_C2, angle1))
            str_CFeC = "C%d-Fe-C%d" %(idx_C1, idx_C2)
            angle_dict[str_CFeC] = angle1
    # find the maximal angle and the corresponding idx
    max_value = max(angle_dict, key=angle_dict.get)
    m1 = angle_dict[max_value]
    str = "Molecule %d %s, %.2f" %(which_mole+1, max_value, angle_dict[max_value])
    angle_dict.pop(max_value, None)
    max_value2 = max(angle_dict, key=angle_dict.get)
    m2 = angle_dict[max_value2]
    str += " next: %s, %.2f" %(max_value2, angle_dict[max_value2])
    str += " diff %.2f exceed thres (30) %d" %(m1 - m2, (m1 - m2) > 30)
    print(str)
    # isolate the index of the axial Cs
    i1 = int(max_value.split("-")[0].split("C")[-1])
    i2 = int(max_value.split("-")[-1].split("C")[-1])
    return [i1, i2]

def get_axial_list(filename = "../qmmm/incav_eq_weak_bath/E0_1e-4_excite_LP_stronger/simu_1.xc.xyz", idx_frame=6):
    data_raw = mda.Universe(filename)
    traj = data_raw.trajectory
    nframes = len(traj)
    print("There are %d frames in %s" %(nframes, filename))
    print("Working under No.%d frame (spacing between frames is 0.25 ps)" %idx_frame)
    print("Time is %.2f ps" %(idx_frame * 0.25))
    print("In the 16 Fe(CO)5 geometries, the maximal C-Fe-C angle is")

    current_pos = traj[idx_frame]._pos.copy()

    xyz_translated = np.zeros((16*11, 3))
    count_change = 0
    axial_idx_lst = []
    for i in range(16):
        nstart = 1911*i
        xyz_feco5 = current_pos[0+nstart:11+nstart]
        xyz_fe = xyz_feco5[0,:]
        xyz_feco5 -= xyz_fe
        # reset the position of each molecule
        xyz_feco5[:,0] += int(i//4) * 8
        xyz_feco5[:,1] += i%4 * 8
        # label if the C is at equatorial or axial locations
        axial_idx = find_axial_Cs(xyz_feco5, which_mole=i)
        axial_idx_lst.append(axial_idx)
    return axial_idx_lst

N_CO_tot = 16 * 5
N_CO_axis = 16 * 2
N_CO_side = 16 * 3
mC_au = 12.0107 * 1822.8885
mO_au = 15.9994 * 1822.8885
Angstrom2AU = 1.8897259885789
au2cminv = 219474.63

# calculate the thermal energy of CO
e_KT = 208.5
e_axis_thermal = e_KT * 2 * 3 * 16
e_side_thermal = e_KT * 3 * 3 * 16

def get_CO_energy_lst(velo):
    # get kinetic energy
    e_lst = []
    for idx in range(1,6):
        v_C = velo[idx,:]
        v_O = velo[idx+5,:]
        e = 0.5 * np.sum(mC_au * v_C**2 + mO_au * v_O**2)
        e_lst.append(e * au2cminv)
    return e_lst

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

def get_qmmm_CO_distr(velo_xyz = "../qmmm/incav_eq_weak_bath/E0_1e-4_excite_LP_stronger/simu_1.out", N=16, window_size=250, t_middle_ps=1.5):
    data = np.loadtxt(velo_xyz)
    nframes = np.shape(data)[0]
    print("There are %d qmmm frames" %nframes)
    e_2d_lst = []
    for idx in range(nframes):
        current_veloc = construct_velo_qmmm(data, idx, N=N)
        e1, e2 = 0.0, 0.0
        es = []
        for idx in range(N):
            current_veloc_local = current_veloc[idx*11:(idx+1)*11,:]
            es += get_CO_energy_lst(current_veloc_local)
        e_2d_lst.append(es)
    e_2d_lst = np.array(e_2d_lst) #- e_KT * 3
    t =  np.array([2.0*i for i in range(nframes)])
    print("kinetic energy of CO bonds have dimension:", np.shape(e_2d_lst))
    print(np.shape(t))
    # The index of the time_middle_ps
    nstep = int(t_middle_ps * 1e3 / 2)
    window_half = int(window_size/2)
    print("average window is %.3f ps" %(window_size*2/1e3))
    e_2d_lst_sliced = e_2d_lst[nstep-window_half:nstep+window_half,:]
    print("sliced kinetic energy of CO bonds have dimension:", np.shape(e_2d_lst_sliced))
    e_1d_lst_sliced = np.mean(e_2d_lst_sliced, axis=0)
    print("sliced  mean kinetic energy of CO bonds have dimension:", np.shape(e_1d_lst_sliced))
    return e_1d_lst_sliced

def filter_CO_distr(e_1d, xyz = "../qmmm/incav_eq_weak_bath/E0_1e-4_excite_LP_stronger/simu_1.xc.xyz", is_axial=True):
    axial_idx_lst = get_axial_list(filename=xyz, idx_frame=6)
    print(axial_idx_lst)
    axial_mapping_lst = []
    if is_axial:
        logic=[True, False]
    else:
        logic=[False, True]
    for i in range(16):
        axial_sub_lst = axial_idx_lst[i]
        for j in range(1, 6):
            if j in axial_sub_lst:
                axial_mapping_lst.append(logic[0])
            else:
                axial_mapping_lst.append(logic[1])
    return axial_mapping_lst


labels=["LP pumping", "UP pumping", "thermal"]
paths = ["../qmmm/incav_eq_weak_bath/E0_1e-4_excite_LP_stronger",
         "../qmmm/incav_eq_weak_bath/E0_1e-4_excite_UP_stronger",
         "../qmmm/incav_eq_weak_bath/E0_1e-4_rerun"]
paths = paths[::-1]
labels = labels[::-1]

axes = clp.initialize(col=3, row=1, width=4.3, height=4.3*0.618*2, fontsize=12,
    sharex=True, sharey=True, LaTeX=True, labelthem=True, labelthemPosition=[-0.1, 1.1])

for i in range(3):
    # preparing three arrays
    e_1d = get_qmmm_CO_distr(velo_xyz=paths[i] + "/simu_1.out")
    axial_mapping_lst = filter_CO_distr(e_1d, xyz=paths[i] + "/simu_1.xc.xyz", is_axial=True)
    equatorial_mapping_lst = filter_CO_distr(e_1d, xyz=paths[i] + "/simu_1.xc.xyz", is_axial=False)
    e_1d_axial = e_1d[axial_mapping_lst]
    e_1d_equatorial = e_1d[equatorial_mapping_lst]

    # calculate energy ratio between two arrays
    e_axial_tot = np.sum(e_1d_axial)
    e_equatorial_tot = np.sum(e_1d_equatorial)
    ratio = e_equatorial_tot / e_axial_tot

    axes[i].hist(e_1d_equatorial, range=(0, 13000), bins=30, density=False, log=False, color = "red", ec="red", alpha=0.5)
    axes[i].hist(e_1d_axial, range=(0, 13000), bins=30, density=False, log=False, color = "0.5", ec="0.5", alpha=0.5)
    axes[i].text(0.15, 0.67, labels[i] + "\nKE$_{\\rm eq}$/KE$_{\\rm ax}$ = %.2f" %ratio, fontsize=12, transform=axes[i].transAxes)

    axes[i].set_xlim(0, 12000)

    axes[i].set_ylabel("CO counts")

#create legend
handles = [Rectangle((0,0),1,1,color=c,ec=c, alpha=0.5) for c in ["red", "0.5"] ]
labels= ["equatorial CO","axial CO"]
axes[0].legend(handles, labels)

axes[-1].set_xlabel("average CO KE [cm$^{-1}$]")
clp.adjust(tight_layout=True, savefile="VCO_dist_LP.pdf")
