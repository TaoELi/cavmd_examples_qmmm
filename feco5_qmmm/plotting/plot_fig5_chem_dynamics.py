import numpy as np
import MDAnalysis as mda
import columnplots as clp
import pandas as pd
import sys

def get_difference(lst1, lst2):
    count_diff = 0
    for idx, item1 in enumerate(lst1):
        item2 = lst2[idx]
        if item1 != item2:
            print("No. %d molecule changes its geometry" %(idx+1))
            count_diff += 1
    if count_diff == 0:
        print("No molecule changes its geometry")
    return count_diff

def C12_at_axial(xyz_feco5, which_mole=0):
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
    #print(str)
    if max_value != "C1-Fe-C2":
        return 1, max_value
    else:
        return 0, max_value

def get_xyz(filename = "../qmmm/incav_eq_weak_bath/E0_1e-4_excite_UP_weaker/simu_1.xc.xyz", idx_frame=0,
            output=None):
    data_raw = mda.Universe(filename)
    traj = data_raw.trajectory
    nframes = len(traj)
    #print("There are %d frames in %s" %(nframes, filename))
    #print("Working under No.%d frame (spacing between frames is 0.25 ps)" %idx_frame)
    print("Time is %.2f ps" %(idx_frame * 0.25))
    #print("In the 16 Fe(CO)5 geometries, the maximal C-Fe-C angle is")

    current_pos = traj[idx_frame]._pos.copy()

    # We now get the 16 Fe(CO)5 geometry at that frame, and also move the center of Fe accordingly
    feco5_str1 = ["Fe", "F", "F", "C", "C", "C", "O", "O", "O", "O", "O"]
    feco5_str2 = ["Fe", "C", "C", "F", "C", "F", "O", "O", "O", "O", "O"]
    feco5_str3 = ["Fe", "C", "C", "C", "F", "F", "O", "O", "O", "O", "O"]
    feco5_str_tot = feco5_str1*9 + feco5_str2 + feco5_str1*2 + feco5_str2*5
    xyz_translated = np.zeros((16*11, 3))
    count_change = 0
    max_value_lst = []
    for i in range(16):
        nstart = 1911*i
        xyz_feco5 = current_pos[0+nstart:11+nstart]
        xyz_fe = xyz_feco5[0,:]
        xyz_feco5 -= xyz_fe
        # reset the position of each molecule
        xyz_feco5[:,0] += int(i//4) * 8
        xyz_feco5[:,1] += i%4 * 8
        # determine the position of the first and second C atoms (labeled as "F")
        count, max_value = C12_at_axial(xyz_feco5, which_mole=i)
        count_change += count
        max_value_lst.append(max_value)
        xyz_translated[0+11*i:11+11*i,:] = xyz_feco5

    #print("%d molecules changed" %count_change)

    xyz = "%d\n\n" %(11*16)
    for i in range(11*16):
        x, y, z = xyz_translated[i,:]
        xyz += "%s %.6f %.6f %.6f\n" %(feco5_str_tot[i], x, y, z)

    if output is not None:
        with open(output, "w") as f:
            f.write(xyz)
    # return the list
    return max_value_lst

def get_qm_xyz(filename = "../qm_solvent_free/incav_eq_weak_bath/E0_1e-4_excite_LP_NVE/simu_1.xc.xyz", idx_frame=0,
            output=None):
    data_raw = mda.Universe(filename)
    traj = data_raw.trajectory
    nframes = len(traj)
    #print("There are %d frames in %s" %(nframes, filename))
    #print("Working under No.%d frame (spacing between frames is 2 fs)" %idx_frame)
    print("Time is %.2f ps" %(idx_frame * 0.002))
    #print("In the 16 Fe(CO)5 geometries, the maximal C-Fe-C angle is")

    current_pos = traj[idx_frame]._pos.copy()

    # We now get the 16 Fe(CO)5 geometry at that frame, and also move the center of Fe accordingly
    feco5_str1 = ["Fe", "F", "F", "C", "C", "C", "O", "O", "O", "O", "O"]
    feco5_str2 = ["Fe", "C", "C", "F", "F", "C", "O", "O", "O", "O", "O"]
    feco5_str_tot = feco5_str1*1 + feco5_str2 + feco5_str1*14
    xyz_translated = np.zeros((16*11, 3))
    count_change = 0
    max_value_lst = []
    for i in range(16):
        nstart = 11*i
        xyz_feco5 = current_pos[0+nstart:11+nstart]
        xyz_fe = xyz_feco5[0,:]
        xyz_feco5 -= xyz_fe
        # reset the position of each molecule
        xyz_feco5[:,0] += int(i//4) * 8
        xyz_feco5[:,1] += i%4 * 8
        # determine the position of the first and second C atoms (labeled as "F")
        count, max_value = C12_at_axial(xyz_feco5)
        count_change += count
        max_value_lst.append(max_value)
        xyz_translated[0+11*i:11+11*i,:] = xyz_feco5

    #print("%d molecules changed" %count_change)

    xyz = "%d\n\n" %(11*16)
    for i in range(11*16):
        x, y, z = xyz_translated[i,:]
        xyz += "%s %.6f %.6f %.6f\n" %(feco5_str_tot[i], x, y, z)

    if output is not None:
        with open(output, "w") as f:
            f.write(xyz)
    # return the list
    return max_value_lst

# plot the number of chemical events as a function of time

UP_accu_lst = [0]
LP_accu_lst = [0]
outcav_accu_lst = [0]
incav_accu_lst = [0]

for idx in range(0, 8, 2):
    # Axial CO for No.idx-1 frame
    lst_UP_p = get_xyz(filename = "../qmmm/incav_eq_weak_bath/E0_1e-4_excite_UP_stronger/simu_1.xc.xyz", idx_frame=idx)
    # Axial CO for No.idx frame
    lst_UP = get_xyz(filename = "../qmmm/incav_eq_weak_bath/E0_1e-4_excite_UP_stronger/simu_1.xc.xyz", idx_frame=idx+2)
    # compare how many reaction events occur between the neighboring frames
    n_UP_diff = get_difference(lst_UP_p, lst_UP)
    # record the accumulated number of events up to this time
    n_UP_p = UP_accu_lst[-1]
    UP_accu_lst.append(n_UP_p + n_UP_diff)

print("\n")
for idx in range(0, 8, 2):
    # Axial CO for No.idx-1 frame
    lst_LP_p = get_xyz(filename = "../qmmm/incav_eq_weak_bath/E0_1e-4_excite_LP_stronger/simu_1.xc.xyz", idx_frame=idx)
    # Axial CO for No.idx frame
    lst_LP = get_xyz(filename = "../qmmm/incav_eq_weak_bath/E0_1e-4_excite_LP_stronger/simu_1.xc.xyz", idx_frame=idx+2)
    # compare how many reaction events occur between the neighboring frames
    n_LP_diff = get_difference(lst_LP_p, lst_LP)
    # record the accumulated number of events up to this time
    n_LP_p = LP_accu_lst[-1]
    LP_accu_lst.append(n_LP_p + n_LP_diff)

print("\n")
for idx in range(0, 8, 2):
    # Axial CO for No.idx-1 frame
    lst_out_p = get_xyz(filename = "../qmmm/incav_eq_weak_bath/outcav_rerun/simu_1.xc.xyz", idx_frame=idx)
    # Axial CO for No.idx frame
    lst_out = get_xyz(filename = "../qmmm/incav_eq_weak_bath/outcav_rerun/simu_1.xc.xyz", idx_frame=idx+2)
    # compare how many reaction events occur between the neighboring frames
    n_outcav_diff = get_difference(lst_out_p, lst_out)
    # record the accumulated number of events up to this time
    n_outcav_p = outcav_accu_lst[-1]
    outcav_accu_lst.append(n_outcav_p + n_outcav_diff)

print("\n")
for idx in range(0, 8, 2):
    # Axial CO for No.idx-1 frame
    lst_in_p = get_xyz(filename = "../qmmm/incav_eq_weak_bath/E0_1e-4_rerun/simu_1.xc.xyz", idx_frame=idx)
    # Axial CO for No.idx frame
    lst_in = get_xyz(filename = "../qmmm/incav_eq_weak_bath/E0_1e-4_rerun/simu_1.xc.xyz", idx_frame=idx+2)
    # compare how many reaction events occur between the neighboring frames
    n_incav_diff = get_difference(lst_in_p, lst_in)
    # record the accumulated number of events up to this time
    n_incav_p = incav_accu_lst[-1]
    incav_accu_lst.append(n_incav_p + n_incav_diff)

UP_accu_lst = np.array(UP_accu_lst)
LP_accu_lst = np.array(LP_accu_lst)
outcav_accu_lst = np.array(outcav_accu_lst)
incav_accu_lst = np.array(incav_accu_lst)

time_lst = np.array([0.5*i for i in range(len(UP_accu_lst))])

xs = [time_lst]*4
ys = [UP_accu_lst, LP_accu_lst, outcav_accu_lst, incav_accu_lst]

colors = ["r--o", "b--s", "k-", "m--"]
labels = ["excite UP", "excite LP", "outside cavity no pulse", "inside cavity no pulse"]

import columnplots as clp
ax = clp.initialize(1, 1, LaTeX=True, width=4.3, fontsize=12)
clp.plotone(xs, ys, ax, colors=colors, labels=labels, xlabel="time [ps]", ylabel="\# barrier crossing events", xlim=[0, 2], ylim=[0,10.5], markersize=10)

# finally, we insert some images
import matplotlib.image as mpimg
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox
data_img = mpimg.imread('feco5_reaction_demo.tiff')
imagebox = OffsetImage(data_img, zoom=0.063)
ab = AnnotationBbox(imagebox, (1.3, 2.2), frameon=False)
ab.set(zorder=-1)
ax.add_artist(ab)

clp.adjust(tight_layout=True, savefile="chem_dynamics_accu.pdf")
