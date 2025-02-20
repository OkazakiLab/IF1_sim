# 4. 120° -> 240° simulation

The simulation files (.pdb, .prmtop, .prmcrd) are stored in ```./4-2_tleap```.

## 4-1. Nucleotide replacement

To see the conformational change of IF1 for dissociation, we continued the simulation by using the snapshot at 100 ns of the γ-stall simulation.  

- αβ1 (IF1-bound, Initialy DP state) should take TP state with ATP.   
- αβ2 (Initialy TP state) should take E state with no nucleotides.  
- αβ3 (Initialy E state) should take DP state with ADP + Pi.  

### ADP+Pi alignment for αβ3 state

We used the βE of the 1h8e structure becuase it has ADP+Pi in the αβE.  
Before alignment, the missing residues of the βE in 1h8e were modeled by MODELLER.  
The modeled 1h8e βE and the simulation structure (β3) were aligned by the Cα atoms of the 9 nucleotide-bound residues; residues 161 162 163 164 189 345 421 424 425.

### Deleting nucleotide from αβ2 state

Just deleting the nucleotide information from the simulation structure.

### Substitution of ADP + Pi into ATP in the αβ1 state

Just deleting Pi and changing the name of 'ADP' to 'ATP'. The following tleap process can add the γ-phosphate, according to the name of 'ATP'. 

## 4-2. tleap

See ```../2_Simulation_setup/2-3_tleap``` for details.

## 4-3. Minimization, Equilibrium run with noH, Cα atoms restraints

See ```../3_0deg_to_120deg_simulation``` for details. Do not forget to 
- add the periodic boundary conditions in the ```min-heat-rnoH.conf``` file
- prepare the restrain files

## 4-4. Equilibrium run with γ angle restrained

Before targeted MD, We run 10 ns γ-stall simulation. 

## 4-5. Targeted MD

To see the complete conformational change of each αβ, 50 ns targeted MD toward αβ2 and αβ3 pairs was performed by NAMD, reducing the root-mean-square error (RMSD) between the current coordinates and the target structure. The target structure for reference was taken from the ground-state structure of bMF1 (PDB: 2jdi): the αβ2 was forced to take the Empty structure, and the αβ3 was forced to take the DP structure. This procedure was not performed on the αβ1 since an appropriate conformational change of the αβ1 was observed during the 0° to 120° γ rotation. 

### Reference file preparation

(1) The missing residues in 2jdi were modeled by MODELLER.  
(2) Each αβ structure in 2jdi was aligned and docked to the simulation structure (after tleap process, ```F1-IF1DP_syn120deg.ff14SB_ionwat.pdb```).    
(3) The resultant (aligned and docked) structure was subject to tleap to add the hydrogen atoms because the original 2jdi structure does not possess them. We do not have to add the ions and waters (see below).  
(4) The number of atoms in the reference file MUST be the same as that of the simulation file ```F1-IF1DP_syn120deg.ff14SB_ionwat.pdb```. The coordinates of nucleotides, ions, and water molecules in the simulation file should be transferred to the reference file. 
(5) Use VMD to set the occupancy value as 1 for the targeted Cα residues (αβ2 and αβ3), and 0 for the non-restrained residues. You should confirm the number of atoms in the reference file again.

## 4-6. Equilibrium run with γ angle restrained (after targeted MD)

The equilibrium run was performed with the γ subunit coordinates restrained for 100 ns.

## 4-7. 120° rotation (120° -> 240°)

See ```../3_0deg_to_120deg_simulation``` for details.

## 4-8. γ-stall at 240° 

After reaching 240° state, 100 ns γ-stall simulation was performed. 