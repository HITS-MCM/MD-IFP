# MD-trajectory_analysis
## A workflow for generation and analysis Protein-Ligand Interaction Fingerprints from MD tajectories


__Authors:__

* Daria Kokh
* Fabian Ormersbach - preprocessing PDB files using Chimera (Process_pdb.py, chimera_hydrogen_mol2.py; test examples revised) 

v.1.0

Daria.Kokh@h-its.org

Heidelberg Institute of Theoretical Studies (HITS, www.h-its.org)

Schloss-Wolfsbrunnenweg 35

69118 Heidelberg, Germany
    

*This open source software code was developed in part in the __Human Brain Project__, funded from the European Union’s Horizon 2020 Framework Programme for Research and Innovation under Specific Grant Agreements  No. 785907 (Human Brain Project  SGA2).*

## __Packages requirements:__
__Python 3.x__

__Python Libraries:__ numpy;    pandas;  matplotlib;  MDAnalysis;  seaborn; RDkit; sklearn;  scipy; ngview 

__Chimera__ - only for the scripts used for preprocessing pdb files (structure protonation and generation of the ligand mol2 file); not requared if protonation and mol2 file are already prepared by a user)
    
__Codes were written on Python 3.x and tested on Python 3.7__

__To configure envirement in anaconda use:__
conda env create -f IFP_trajectory.yml



## Scripts:


__Trajectories.py__  - functions for building a trajectory object for reading and analysis standard MD and RAMD trajectories and computation of relative residence times

__IFP_generation.py__  -  functions for generation of IFPs

__Membrane_analysis.py__ - functions for analysis of membrane-protein systems 

__Clustering.py__   - functions for analysis of trajectories using IFP data   (is still under developments)

__Process_pdb.py__   - preprocessing PDB files (splitting into ligand and protein files)

__chimera_hydrogen_mol2.py__  - generation ligand mole2 file

__Membrane_analysis.py__ - computation of the membrane surface area per lipid and membrane/protein/water atom density distribution  


## Data employed in test examples: 
https://zenodo.org/record/3755337#.XrF-iGgzaUk
       
## Test example as jupyter notebooks :

I. __IFP_generation_examples_PDB.ipynb:__

Protein-Ligand interaction fingerprint (IFP) computations (only function of IFP_generation.py are used) for:
   1. a single structure prepared for MD simulations (HSP90; PDB ID 6EI5, dcd format)
   2. a trajectory (for selected frames; dcd format)
   3. a PDB structure

II. __IFP_generation_examples_TRAJ.ipynb:__ 

generation and analysis of  IFPs for plain MD simulations and RAMD trajectories for Muscarinic Receptor M2 in a membrane
in this example Trajectories.py is used for pre-processing trajectories and  IFP_generation.py for computing IFPs
   1. Computing IFPs for system equilibration trajectory (dcd format)
   3. Computing IFPs for ligand dissociation trajectories (dcd format)

III. __IFP_generation_examples_Analysis.ipynb:__ 

This example shows how RAMD dissociation trajectories can be analyzed usinf pre-generated IFP database 
![HSP90](/images/cluster-traj.png)
*This plot illustrates ligand dissociation pathways in a graph representation*
   
IV. __membrane_analysis_example.ipynb__

Exploring behavior of a membrane-protein system in MD simulations

