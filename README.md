# MD-trajectory_analysis
## A workflow for generation and analysis Protein-Ligand Interaction Fingerprints from MD tajectories


__Authors:__

Daria Kokh

Fabian Ormersbach (chimera_hydrogen_mol2.py; test examples revised) 

v.1.0
Daria.Kokh@h-its.org
Heidelberg Institute of Theoretical Studies (HITS, www.h-its.org)
Schloss-Wolfsbrunnenweg 35
69118 Heidelberg, Germany
    

*This open source software code was developed in part in the Human Brain Project, funded from the European Union’s Horizon 2020 Framework Programme for Research and Innovation under Specific Grant Agreements  No. 785907 (Human Brain Project  SGA2).*

## __Packages requirements:__
    Python 3.x
    numpy;    pandas;  matplotlib;  MDAnalysis;  RDkit;   scipy; ngview
    for 
    
__Codes were written on Python 3.x and tested on Python 3.7__
__To configure envirement in anaconda use__
conda env create -f IFP_trajectory.yml



## Scripts:

__Clustering.py__   - functions for analysis of trajectories using IFP data   (is still under developments)

__IFP_generation.py__  -  functions for generation of IFPs
__Membrane_analysis.py__ - functions for analysis of membrane-protein systems 
__Trajectories.py__  - functions for building a trajectory object for reading and analysis standard MD and RAMD trajectories and computation of relative residence times

## Data employed in test examples: 
https://zenodo.org/record/3755337#.XrF-iGgzaUk
       
## Example jupyter notebooks :

1. __IFP_generation_examples_PDB.ipynb:__

an examples of Protein-Ligand IFP computations for
  ** a single structure prepared for MD simulations (HSP90; PDB ID 6EI5, dcd format)
  ** a trajectory (for selected frames; dcd format)
  ** a PDB structure

2. __IFP_generation_example_TRAJ.ipynb: example of generation and analysis of Ligand-Protein IFPs for RAMD simulations (dcd trajectories):__
  ** Computing interaction fingerprints (IFP) for
  ** system equilibration trajectory
  ** ligand dissociation trajectories
