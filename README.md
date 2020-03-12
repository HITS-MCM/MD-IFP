# MD-trajectory_analysis
Python scripts for analysis of Molecular Dynamics trajectories (based on the MDAnalysis and RDKit libraries).
    v.1.0
    Copyright (c) 2019
    Released under the GNU Public Licence, v2 or any higher version


This open source software code was developed in part in the Human Brain Project, funded from the European Union’s Horizon 2020 Framework Programme for Research and Innovation under Specific Grant Agreements  No. 785907 (Human Brain Project  SGA2).


Packages required:
    numpy;    pandas;  matplotlib;  MDAnalysis;  RDkit;   scipy
    code is written on Python 3.x and tested on Python 3.7
   
    to configure envirement in anaconda use
    conda env create -f IFP_trajectory.yml


Scripts:

1. Trajectories.py 	- analysis of protein-ligand dissociation pathways
2. Membrane_analysis.py - analysis of the membrane properties
3. IFP_generation.py  	- generation of Protein-Ligand Interaction Fingerprints  (PL IFP)


Examples:
1. Example jupyter notebooks : 

  (i)IFP_generation_examples_PDB.ipynb
   Test examples of PL IFP computations
   1. Computing interaction fingerprints (IFP) for
     -- a single structure prepared for MD simulations (HSP90; PDB ID 6EI5, dcd format)
     -- a trajectory (for selected frames; dcd format)
     -- PDB structure
   2. Visualizing protein residues that are involved in protein-ligand interactions, including water-bridges

  (ii)IFP_generation_example_TRAJ.ipynb
   Generation and analysis of Ligand-Protein IFPs for RAMD simulations (dcd trajectories)
   1. Computing interaction fingerprints (IFP) for 
      system equilibration trajectory 
      ligand dissociation trajectories
   2. Visualizion of 
      protein residues that are involved in protein-ligand interactions, including water-bridges
      ligand dissociation   

2. Data:
   	2YKI  
	4MQT  
	6EI5 

