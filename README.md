# MD-IFP: MD trajectory analysis using protein-ligand Interaction Fingerprints
## A Python Workflow for the Generation and Analysis of Protein-Ligand Interaction Fingerprints from Molecular Dynamics trajectories
## v.1.0
## 03.05.2020

__Authors and Contributors:__

* Daria Kokh
* Fabian Ormersbach - preprocessing PDB files using Chimera (Process_pdb.py, chimera_hydrogen_mol2.py; test examples revised) 


Daria.Kokh@h-its.org

Heidelberg Institute of Theoretical Studies (HITS, www.h-its.org)

Schloss-Wolfsbrunnenweg 35

69118 Heidelberg, Germany
    

*This open source software code was developed in part in the __Human Brain Project__, funded from the European Union’s Horizon 2020 Framework Programme for Research and Innovation under Specific Grant Agreements  No. 785907 (Human Brain Project  SGA2).*

## __Packages requirements:__
__Python 3.x__

__Python Libraries:__ numpy;    pandas;  matplotlib;  MDAnalysis;  seaborn; RDkit; sklearn;  scipy; ngview 

__Chimera__ - only for the scripts used for preprocessing pdb files (structure protonation and generation of the ligand mol2 file); not required if protonation and mol2 file are already prepared by a user)
    
__Codes were written on Python 3.x and tested on Python 3.7__

__To configure environment in anaconda use:__
conda env create -f IFP_trajectory.yml



## Scripts:


### __Trajectories.py__  - functions for building a trajectory object for reading and analysis of standard MD and RAMD trajectories and computation of relative residence times

__IFP_generation.py__  -  functions for generation of IFPs

__Membrane_analysis.py__ - functions for analysis of membrane-protein systems 

__Clustering.py__   - functions for analysis of trajectories using IFP data   (is still under developments)

__Process_pdb.py__   - preprocessing PDB files (splitting into ligand and protein files)

__chimera_hydrogen_mol2.py__  - generation of ligand mol2 file

__Membrane_analysis.py__ - computation of the membrane surface area per lipid and membrane/protein/water atom density distribution  

       
## Application examples (folder Examples):

   1. Generation of the IFP databease for a single MD trajectory of a protein-ligand complex
   2. Analysis and visualization of a set of IFP databases for different ligand 

## Test Examples as Python Jupyter Notebooks :

### Data employed in test examples 
   can be downloaded from  https://zenodo.org/record/3755337#.XrF-iGgzaUk

### I. __IFP_generation_examples_PDB.ipynb:__

Protein-Ligand Interaction Fingerprint (IFP) computations (only functions of IFP_generation.py are used) for:
   1. a single structure prepared for MD simulations (HSP90; PDB ID 6EI5, dcd format)
   2. a trajectory (for selected frames; dcd format)
   3. a PDB structure


### II. __IFP_generation_examples_TRAJ.ipynb:__ 

Generation and analysis of IFPs for conventional MD simulations and for RAMD trajectories for Muscarinic Receptor M2 in a membrane.  In this example,  Trajectories.py is used for pre-processing trajectories and IFP_generation.py is used for computing IFPs
   1. Computing IFPs for a single equilibration trajectory (dcd format)
   3. Computing IFPs for a set of trajectories: system equilibration and ligand dissociation (RAMD) trajectories (dcd format)
![HSP90](/images/ifp_RAMD_4MQT.png)
*Illustration of PL IFP variation in one of the dissociation trajectories of iperoxo bound to muscarinic receptor M2 .*


### III. __IFP_generation_examples_Analysis.ipynb:__ 

This example shows how RAMD dissociation trajectories can be analyzed using pre-generated IFP databases 
![HSP90](/images/cluster-traj.png)

*This plot illustrates ligand dissociation pathways in a graph representation derived from clustering ligand trajectories in IFP space and plotting them with respect to RMSD from the initial bound position.*
   
   
### IV. __membrane_analysis_example.ipynb__

Exploring the behavior of a membrane-protein system in MD simulations. The membrane surface is assumed to be in the X/Y plane 
![HSP90](/images/Membrane.png)
*Illustratration of the analysis of the muscarinic receptor M2 GPCR embedded in a mixed membrane with 50% cholesterol content*
