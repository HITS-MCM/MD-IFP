# MD-trajectory_analysis
## A workflow for generation and analysis Protein-Ligand Interaction Fingerprints from MD tajectories


__Authors:__

Daria Kokh

Fabian Ormersbach preprocessing PDB files using Chimera (Process_pdb.py, chimera_hydrogen_mol2.py; test examples revised) 

v.1.0
Daria.Kokh@h-its.org
Heidelberg Institute of Theoretical Studies (HITS, www.h-its.org)
Schloss-Wolfsbrunnenweg 35
69118 Heidelberg, Germany
    

*This open source software code was developed in part in the Human Brain Project, funded from the European Union’s Horizon 2020 Framework Programme for Research and Innovation under Specific Grant Agreements  No. 785907 (Human Brain Project  SGA2).*

## __Packages requirements:__
    Python 3.x
    numpy;    pandas;  matplotlib;  MDAnalysis;  RDkit;   scipy; ngview
    Chimera - only for the scripts used for preprocessing pdb files (structure protonation and generation of the ligand mol2 file); not requared if protonation and mol2 file are already prepared by a user)
    
__Codes were written on Python 3.x and tested on Python 3.7__

__To configure envirement in anaconda use__
conda env create -f IFP_trajectory.yml



## Scripts:


__Trajectories.py__  - functions for building a trajectory object for reading and analysis standard MD and RAMD trajectories and computation of relative residence times

__IFP_generation.py__  -  functions for generation of IFPs

__Membrane_analysis.py__ - functions for analysis of membrane-protein systems 

__Clustering.py__   - functions for analysis of trajectories using IFP data   (is still under developments)



## Data employed in test examples: 
https://zenodo.org/record/3755337#.XrF-iGgzaUk
       
## Test example as jupyter notebooks :

I. __IFP_generation_examples_PDB.ipynb:__

an examples of Protein-Ligand IFP computations for
   1. a single structure prepared for MD simulations (HSP90; PDB ID 6EI5, dcd format)
   2. a trajectory (for selected frames; dcd format)
   3. a PDB structure

II. __IFP_generation_examples_TRAJ.ipynb:__ 

example of generation and analysis of Ligand-Protein IFPs for RAMD simulations (dcd trajectories)

   1. Computing interaction fingerprints (IFP) for
   2. system equilibration trajectory
   3. ligand dissociation trajectories


III. __IFP_generation_examples_Analysis.ipynb:__ 

example of RAMD dissociation trajectories analysis using pre-generated IFP database 

   1. Example of HSP90
