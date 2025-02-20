from audioop import rms
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
import os

os.chdir('/Users/ryohei/desktop')

sim = mda.Universe('F1-IF1DP_syn120deg_stall_ATPchainD.pdb')
b1h8e = mda.Universe('1h8e_E_ADPPi_rename.pdb')

# Chain D (Beta DP): ADP+Pi -> ATP
#------------------------------------------------------------------------------------#
sim_E = sim.select_atoms('(resid 2084 2085 2086 2087 2112 2268 2344 2347 2348) and name CA') #resid 1932 to 2049 or resid 2052 to 2387
sim_E0 = sim_E.positions - sim_E.atoms.center_of_mass()

b1h8e_E = b1h8e.select_atoms('protein and (resid 161 162 163 164 189 345 421 424 425) and name CA')
b1h8e_E0 = b1h8e_E.positions - b1h8e_E.atoms.center_of_mass()
print(len(sim_E0), len(b1h8e_E0))

R, rmsd = align.rotation_matrix(b1h8e_E0, sim_E0)
print(R, rmsd)

b1h8e.atoms.translate(-b1h8e_E.select_atoms('all').center_of_mass())
b1h8e.atoms.rotate(R)
b1h8e.atoms.translate(sim_E.select_atoms('all').center_of_mass())
b1h8e.atoms.write("1h8e_IF1DP_mobile.pdb")

#------------------------------------------------------------------------------------#

#pdb_r = 'F1-noIF1_syn120deg_stall_ATPchainD.pdb'
#b1h8e_r = '1h8e_mobile_D.pdb'
#pdb_w = 'F1-noIF1_syn120deg_stall_nucletide_rep.pdb'

# ValueError: invalid literal for int() with base 10: '     '
# と言うエラーが出るので、1h8eからのヌクレオチドの置換は手動で行う