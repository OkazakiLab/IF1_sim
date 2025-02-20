import numpy as np
import matplotlib.pyplot as plt
import math
import csv
import re
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = 18
#plt.rcParams.update({'font.size': 16})

b1bmf = './bMF1/1bmf_center_rotz_roty_rotz.pdb'
b1cow = './bMF1/1cow_center_rotz_roty_rotz.pdb'
b1e1q = './bMF1/1e1q_center_rotz_roty_rotz.pdb'
b1e1r = './bMF1/1e1r_center_rotz_roty_rotz.pdb'
b1e79 = './bMF1/1e79_center_rotz_roty_rotz.pdb'
b1efr = './bMF1/1efr_center_rotz_roty_rotz.pdb'
b1h8e = './bMF1/1h8e_center_rotz_roty_rotz.pdb'
b1h8h = './bMF1/1h8h_center_rotz_roty_rotz.pdb'
b1nbm = './bMF1/1nbm_center_rotz_roty_rotz.pdb'
b1ohh = './bMF1/1ohh_center_rotz_roty_rotz.pdb'
b1w0j = './bMF1/1w0j_center_rotz_roty_rotz.pdb'
b1w0k = './bMF1/1w0k_center_rotz_roty_rotz.pdb'
b2ck3 = './bMF1/2ck3_center_rotz_roty_rotz.pdb'
b2jdi = './bMF1/2jdi_center_rotz_roty_rotz.pdb'
b2jiz_1 = './bMF1/2jiz_1_center_rotz_roty_rotz.pdb'
b2jiz_2 = './bMF1/2jiz_2_center_rotz_roty_rotz.pdb'
b2jj1_1 = './bMF1/2jj1_1_center_rotz_roty_rotz.pdb'
b2jj1_2 = './bMF1/2jj1_2_center_rotz_roty_rotz.pdb'
b2jj2_1 = './bMF1/2jj2_1_center_rotz_roty_rotz.pdb'
b2jj2_2 = './bMF1/2jj2_2_center_rotz_roty_rotz.pdb'
b2v7q = './bMF1/2v7q_center_rotz_roty_rotz.pdb'
b4asu = './bMF1/4asu_center_rotz_roty_rotz.pdb'
b4tsf = './bMF1/4tsf_center_rotz_roty_rotz.pdb'
b4tt3 = './bMF1/4tt3_center_rotz_roty_rotz.pdb'
b4yxw = './bMF1/4yxw_center_rotz_roty_rotz.pdb'
b4z1m = './bMF1/4z1m_center_rotz_roty_rotz.pdb'

pdblist = [b1bmf, b1cow, b1e1q, b1e1r, b1e79, b1efr, b1h8e, b1h8h, b1nbm, b1ohh, 
           b1w0j, b1w0k, b2ck3, b2jdi, b2jiz_1, b2jiz_2, b2jj1_1, b2jj1_2, b2jj2_1, b2jj2_2, 
           b2v7q, b4asu, b4tsf, b4tt3, b4yxw, b4z1m]
             
ref = mda.Universe(b1bmf).select_atoms('segid G and name CA and resid 1:30, 221:270').positions
L = len(pdblist)

pdb_names = []
thetas = []

with open('Summary.csv','a',newline='') as f:
     writer = csv.writer(f)
     writer.writerow(['F1_species', 'PDB', 'rmsd', 'a(2,1)', 'cosine', 'theta', 'Matrix'])

for l in range(L):
    str = pdblist[l]
    if str.startswith('./b'):
        F1_species = 'bMF1'
        mol = mda.Universe(str).select_atoms('segid G and name CA and resid 1:30, 221:270').positions
        if str.endswith('_1_center_rotz_roty_rotz.pdb') or str.endswith('_2_center_rotz_roty_rotz.pdb'):
            pdb_name = str[7:13]
        else:
            pdb_name = str[7:11]

    R, rmsd = align.rotation_matrix(mol, ref)
    print('rmsd =', rmsd)
    print(R)
    
    ref_vector = [1.0, 0.0, 0.0]
    R1 = np.dot(R, ref_vector)
    vector = [R1[0], R1[1], 0]
    
    vec_norm = np.linalg.norm(vector)
    ref_norm = np.linalg.norm(ref_vector)
    cosine = np.dot(vector, ref_vector)/(vec_norm * ref_norm)
    
    if R1[1] < 0.0:
        theta = np.degrees(np.arccos(cosine))
    else:
        theta = -np.degrees(np.arccos(cosine))
    print('cosine=', cosine)    
    print('theta=',  theta)
    
    thetas.append(theta)
    pdb_names.append(pdb_name)

plt.figure(figsize=(9, 6))
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)

plt.scatter(pdb_names, thetas, color='gray', s=90)  
plt.scatter(pdb_names[9], thetas[9], color='orangered', s=90) 
plt.scatter(pdb_names[20], thetas[20], color='orangered', s=90) 
plt.scatter(pdb_names[22], thetas[22], color='orangered', s=90) 
plt.scatter(pdb_names[23], thetas[23], color='orangered', s=90) 
plt.scatter(pdb_names[25], thetas[25], color='orangered', s=90) 

plt.xlabel('PDB Name')
plt.ylabel('Rotary angle (deg)')
plt.xticks(rotation=90) 
plt.tick_params(direction='in')
plt.tight_layout()

plt.savefig("g-angle_bMF1.png")
plt.show()

with open('g-angle_bMF1.txt', 'w') as txt_file:
    txt_file.write("pdb_name\ttheta\n")
    for pdb_name, theta in zip(pdb_names, thetas):
        txt_file.write(f"{pdb_name}\t{theta:.3f}\n")
