from __future__ import print_function
import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = "17"
plt.tick_params(direction="in", length=6, width=1)
plt.gca().spines["top"].set_visible(False)
plt.gca().spines["right"].set_visible(False)

b1bmf_E = './pdb/1bmf_E_selpca.pdb'
b1bmf_TP = './pdb/1bmf_TP_selpca.pdb'
b1bmf_DP = './pdb/1bmf_DP_selpca.pdb'
b1cow_E = './pdb/1cow_E_selpca.pdb'
b1cow_TP = './pdb/1cow_TP_selpca.pdb'
b1cow_DP = './pdb/1cow_DP_selpca.pdb'
b1e1q_E = './pdb/1e1q_E_selpca.pdb'
b1e1q_TP = './pdb/1e1q_TP_selpca.pdb'
b1e1q_DP = './pdb/1e1q_DP_selpca.pdb'
b1e1r_E = './pdb/1e1r_E_selpca.pdb'
b1e1r_TP = './pdb/1e1r_TP_selpca.pdb'
b1e1r_DP = './pdb/1e1r_DP_selpca.pdb'
b1e79_E = './pdb/1e79_E_selpca.pdb'
b1e79_TP = './pdb/1e79_TP_selpca.pdb'
b1e79_DP = './pdb/1e79_DP_selpca.pdb'
b1efr_E = './pdb/1efr_E_selpca.pdb'
b1efr_TP = './pdb/1efr_TP_selpca.pdb'
b1efr_DP = './pdb/1efr_DP_selpca.pdb'
b1h8e_E = './pdb/1h8e_E_selpca.pdb'
b1h8e_TP = './pdb/1h8e_TP_selpca.pdb'
b1h8e_DP = './pdb/1h8e_DP_selpca.pdb'
b1h8h_E = './pdb/1h8h_E_selpca.pdb'
b1h8h_TP = './pdb/1h8h_TP_selpca.pdb'
b1h8h_DP = './pdb/1h8h_DP_selpca.pdb'
b1nbm_E = './pdb/1nbm_E_selpca.pdb'
b1nbm_TP = './pdb/1nbm_TP_selpca.pdb'
b1nbm_DP = './pdb/1nbm_DP_selpca.pdb'
b1ohh_E = './pdb/1ohh_E_selpca.pdb'
b1ohh_TP = './pdb/1ohh_TP_selpca.pdb'
b1ohh_DP = './pdb/1ohh_DP_selpca.pdb'
b1w0j_E = './pdb/1w0j_E_selpca.pdb'
b1w0j_TP = './pdb/1w0j_TP_selpca.pdb'
b1w0j_DP = './pdb/1w0j_DP_selpca.pdb'
b1w0k_E = './pdb/1w0k_E_selpca.pdb'
b1w0k_TP = './pdb/1w0k_TP_selpca.pdb'
b1w0k_DP = './pdb/1w0k_DP_selpca.pdb'
b2ck3_E = './pdb/2ck3_E_selpca.pdb'
b2ck3_TP = './pdb/2ck3_TP_selpca.pdb'
b2ck3_DP = './pdb/2ck3_DP_selpca.pdb'
b2jdi_E = './pdb/2jdi_E_selpca.pdb'
b2jdi_TP = './pdb/2jdi_TP_selpca.pdb'
b2jdi_DP = './pdb/2jdi_DP_selpca.pdb'
b2jiz_1_E = './pdb/2jiz_1_E_selpca.pdb'
b2jiz_1_TP = './pdb/2jiz_1_TP_selpca.pdb'
b2jiz_1_DP = './pdb/2jiz_1_DP_selpca.pdb'
b2jiz_2_E = './pdb/2jiz_2_E_selpca.pdb'
b2jiz_2_TP = './pdb/2jiz_2_TP_selpca.pdb'
b2jiz_2_DP = './pdb/2jiz_2_DP_selpca.pdb'
b2jj1_1_E = './pdb/2jj1_1_E_selpca.pdb'
b2jj1_1_TP = './pdb/2jj1_1_TP_selpca.pdb'
b2jj1_1_DP = './pdb/2jj1_1_DP_selpca.pdb'
b2jj1_2_E = './pdb/2jj1_2_E_selpca.pdb'
b2jj1_2_TP = './pdb/2jj1_2_TP_selpca.pdb'
b2jj1_2_DP = './pdb/2jj1_2_DP_selpca.pdb'
b2jj2_1_E = './pdb/2jj2_1_E_selpca.pdb'
b2jj2_1_TP = './pdb/2jj2_1_TP_selpca.pdb'
b2jj2_1_DP = './pdb/2jj2_1_DP_selpca.pdb'
b2jj2_2_E = './pdb/2jj2_2_E_selpca.pdb'
b2jj2_2_TP = './pdb/2jj2_2_TP_selpca.pdb'
b2jj2_2_DP = './pdb/2jj2_2_DP_selpca.pdb'
b2v7q_E = './pdb/2v7q_E_selpca.pdb'
b2v7q_TP = './pdb/2v7q_TP_selpca.pdb'
b2v7q_DP = './pdb/2v7q_DP_selpca.pdb'
b4asu_E = './pdb/4asu_E_selpca.pdb'
b4asu_TP = './pdb/4asu_TP_selpca.pdb'
b4asu_DP = './pdb/4asu_DP_selpca.pdb'
b4tsf_E = './pdb/4tsf_E_selpca.pdb'
b4tsf_TP = './pdb/4tsf_TP_selpca.pdb'
b4tsf_DP = './pdb/4tsf_DP_selpca.pdb'
b4tt3_E = './pdb/4tt3_E_selpca.pdb'
b4tt3_TP = './pdb/4tt3_TP_selpca.pdb'
b4tt3_DP = './pdb/4tt3_DP_selpca.pdb'
b4yxw_E = './pdb/4yxw_E_selpca.pdb'
b4yxw_TP = './pdb/4yxw_TP_selpca.pdb'
b4yxw_DP = './pdb/4yxw_DP_selpca.pdb'
b4z1m_E = './pdb/4z1m_E_selpca.pdb'
b4z1m_TP = './pdb/4z1m_TP_selpca.pdb'
b4z1m_DP = './pdb/4z1m_DP_selpca.pdb'

structures = [b1bmf_E, b1bmf_TP, b1bmf_DP, b1cow_E, b1cow_TP, b1cow_DP, b1e1q_E, b1e1q_TP, b1e1q_DP,
              b1e1r_E, b1e1r_TP, b1e1r_DP, b1e79_E, b1e79_TP, b1e79_DP, b1efr_E, b1efr_TP, b1efr_DP, 
              b1h8e_E, b1h8e_TP, b1h8e_DP, b1h8h_E, b1h8h_TP, b1h8h_DP, b1nbm_E, b1nbm_TP, b1nbm_DP, 
              b1ohh_E, b1ohh_TP, b1ohh_DP, b1w0j_E, b1w0j_TP, b1w0j_DP, b1w0k_E, b1w0k_TP, b1w0k_DP,
              b2ck3_E, b2ck3_TP, b2ck3_DP, b2jdi_E, b2jdi_TP, b2jdi_DP, b2jiz_1_E, b2jiz_1_TP, b2jiz_1_DP, 
              b2jiz_2_E, b2jiz_2_TP, b2jiz_2_DP, b2jj1_1_E, b2jj1_1_TP, b2jj1_1_DP, b2jj1_2_E, b2jj1_2_TP, b2jj1_2_DP, 
              b2jj2_1_E, b2jj2_1_TP, b2jj2_1_DP, b2jj2_2_E, b2jj2_2_TP, b2jj2_2_DP, b2v7q_E, b2v7q_TP, b2v7q_DP,
              b4asu_E, b4asu_TP, b4asu_DP, b4tsf_E, b4tsf_TP, b4tsf_DP, b4tt3_E, b4tt3_TP, b4tt3_DP,
              b4yxw_E, b4yxw_TP, b4yxw_DP, b4z1m_E, b4z1m_TP, b4z1m_DP] 

pdb_files = md.load(structures)
for i in range(5):
    if i == 0:
        pdb_files.superpose(pdb_files, 0)
    else:
        avs = md.load('average_structure.pdb')
        pdb_files.superpose(avs)
    average = pdb_files.xyz.mean(axis=0)
    topology = pdb_files.top
    average_structure = md.Trajectory(average, topology)
    average_structure.save('average_structure.pdb')

pdb_files.save_pdb('aligned.pdb')
traj = md.load('aligned.pdb')
coords = traj.xyz.reshape(traj.xyz.shape[0], -1)
#print(coords)
cov_matrix = np.cov(coords.T)
#print(cov_matrix)

evals, evecs = np.linalg.eigh(cov_matrix)
#print(evals[-10:]) # Print the 10 most significant eigenvalues.
#print(evecs[:, -1].shape)
evals_sum = 0
for i in range(10):
    evals_sum += evals[-i]

color_list = ['red', 'blue', 'green', 'red', 'blue', 'green', 'red', 'blue', 'green', 'red', 'blue', 'green',
              'red', 'blue', 'green', 'red', 'blue', 'green', 'red', 'blue', 'green', 'red', 'blue', 'green',
              'red', 'blue', 'green', 'red', 'blue', 'green', 'red', 'blue', 'green', 'red', 'blue', 'green',
              'red', 'blue', 'green', 'red', 'blue', 'green', 'red', 'blue', 'green', 'red', 'blue', 'green',
              'red', 'blue', 'green', 'red', 'blue', 'green', 'red', 'blue', 'green', 'red', 'blue', 'green',
              'red', 'blue', 'green', 'red', 'blue', 'green', 'red', 'blue', 'green', 'red', 'blue', 'green',
              'red', 'blue', 'green', 'red', 'blue', 'green']
labels = ['1bmf_E', '1bmf_TP', '1bmf_DP', '1cow_E', '1cow_TP', '1cow_DP', '1e1q_E', '1e1q_TP', '1e1q_DP',
          '1e1r_E', '1e1r_TP', '1e1r_DP', '1e79_E', '1e79_TP', '1e79_DP', '1efr_E', '1efr_TP', '1efr_DP', 
          '1h8e_E', '1h8e_TP', '1h8e_DP', '1h8h_E', '1h8h_TP', '1h8h_DP', '1nbm_E', '1nbm_TP', '1nbm_DP', 
          '1ohh_E', '1ohh_TP', '1ohh_DP', '1w0j_E', '1w0j_TP', '1w0j_DP', '1w0k_E', '1w0k_TP', '1w0k_DP',
          '2ck3_E', '2ck3_TP', '2ck3_DP', '2jdi_E', '2jdi_TP', '2jdi_DP', '2jiz_1_E', '2jiz_1_TP', '2jiz_1_DP',
          '2jiz_2_E', '2jiz_2_TP', '2jiz_2_DP', '2jj1_1_E', '2jj1_1_TP', '2jj1_1_DP', '2jj1_2_E', '2jj1_2_TP', '2jj1_2_DP', 
          '2jj2_1_E', '2jj2_1_TP', '2jj2_1_DP', '2jj2_2_E', '2jj2_2_TP', '2jj2_2_DP', '2v7q_E', '2v7q_TP', '2v7q_DP',
          '4asu_E', '4asu_TP', '4asu_DP', '4tsf_E', '4tsf_TP', '4tsf_DP', '4tt3_E', '4tt3_TP', '4tt3_DP',
          '4yxw_E', '4yxw_TP', '4yxw_DP', '4z1m_E', '4z1m_TP', '4z1m_DP'] 

avs = coords.mean(axis=0)# Compute the average geometry
print(avs.reshape(1, -1))

# Plot projections on the two-most significant principal modes.
proj1 = np.dot(coords - avs.reshape(1, -1), evecs[:, -1])
proj2 = np.dot(coords - avs.reshape(1, -1), evecs[:, -2])
PC1 = evals[-1]/evals_sum
PC2 = evals[-2]/evals_sum

print('PC1=', '{:.3f}'.format(PC1), 'PC2=', '{:.3f}'.format(PC2))
plt.scatter(proj1, proj2, c=color_list, s=50)
plt.xlabel("PC1 ({:.1%})".format(PC1))
plt.ylabel("PC2 ({:.1%})".format(PC2))
plt.xticks([-10, -8, -6, -4, -2, 0, 2, 4, 6]) 
plt.yticks([-4, -2, 0, 2, 4])
plt.savefig("PCA.pdf")
#plt.show()

output_data = np.column_stack((proj1, proj2))
np.savetxt('pc1_pc2.txt', output_data, fmt='%.3f', header='PC1 PC2', comments='')

np.savetxt("example_eigenvector1.txt", evecs[:, -1])
np.savetxt("example_eigenvector2.txt", evecs[:, -2])