# MD-IFP: MD trajectory analysis using protein-ligand Interaction Fingerprints
## A Python Workflow for the Generation and Analysis of Protein-Ligand Interaction Fingerprints from Molecular Dynamics trajectories
### v.1.1
### 19.06.2021
   
## Associated data: 
https://zenodo.org/record/3981155#.XzQEUCgzaUk

## Publications describing the IFP analysis. Please cited the following paper : 
   D. B. Kokh, B. Doser, S. Richter, F. Ormersbach, X. Cheng , R.C. Wade  "A Workflow for Exploring Ligand Dissociation    from a Macromolecule: Efficient Random Acceleration Molecular Dynamics Simulation and Interaction Fingerprints Analysis of Ligand Trajectories" J. Chem Phys.(2020) 158  125102  doi: 10.1063/5.0019088; 
   
## Publications of application examples: 
1. IFP analysis of dissociation trajectories for 3 compounds of HSP90  reported in the paper 
  
   D. B. Kokh, B. Doser, S. Richter, F. Ormersbach, X. Cheng , R.C. Wade  "A Workflow for Exploring Ligand Dissociation    from a Macromolecule: Efficient Random Acceleration Molecular Dynamics Simulation and Interaction Fingerprints Analysis of Ligand         Trajectories" J. Chem Phys.(2020) 153  125102  doi: 10.1063/5.0019088; https://arxiv.org/abs/2006.11066
   
   Results are implemented in  __IFP_generation_examples_Analysis.ipynb__ 
  
2.  Small compound unbinding from T4 lysozyme mutants

    A Nunes-Alves, DB Kokh, RC Wade  "Ligand unbinding mechanisms and kinetics for T4 lysozyme mutants from τRAMD simulations", Current Research in Structural Biology 3, 106-111
    https://doi.org/10.1016/j.crstbi.2021.04.001

3. Application to two GPCR targets (embedded in a membrane):

   D. B. Kokh, R.C. Wade "G-Protein Coupled Receptor-Ligand Dissociation Rates and Mechanisms from τRAMD Simulations",    doi: https://doi.org/10.1101/2021.06.20.449151

   Associated scripts and data can be downloaded here: https://zenodo.org/record/5001884#.YM-rRmgzYuU
   

## Tutorials: 
1. Youtube lecture/tutorial for 2020 MolSSI School on Open Source Software in Rare Event Path Sampling Strategies: "tauRAMD workflow: fast estimation of protein-ligand residence times with insights into dissociation mechanisms" : https://www.youtube.com/watch?v=kCUyQtoo4cE&feature=youtu.be


##  Authors and Contributors:

* Daria Kokh
* Fabian Ormersbach - preprocessing PDB files using Chimera (Process_pdb.py, chimera_hydrogen_mol2.py; test examples revised) 

Daria.Kokh@h-its.org

Heidelberg Institute of Theoretical Studies (HITS, www.h-its.org)

Schloss-Wolfsbrunnenweg 35

69118 Heidelberg, Germany
    
*This open source software code was developed in part in the __Human Brain Project__, funded from the European Union’s Horizon 2020 Framework Programme for Research and Innovation under Specific Grant Agreements  No. 785907 (Human Brain Project  SGA2).*

## __Packages requirements:__
__Python 3.x__ 

__Python Libraries:__ 
   1. numpy;    pandas;  matplotlib;  seaborn; sklearn;  scipy; 
   2. __RDkit__ 
   3. __ngview__ - used for visualization (installation of ngview can be tricky, the following way may work: after installation of the Python envirenment - conda install -c conda-forge nglview=2.7.1 and then jupyter-nbextension enable nglview --py --sys-prefix). If you don't need visualization, you can skip this, but JN must be edited accordingly
   4. __MDAnalysis Version: 0.20.1__ (Important: an old module for H-bond analysis is currently used, it will be removed in version 2.0 )

__Chimera__ - only for the scripts used for preprocessing pdb files (structure protonation and generation of the ligand mol2 file); not required if protonation and mol2 file are already prepared by a user)
    
__Codes were written on Python 3.x and tested on Python 3.7__

__To configure environment in anaconda use:__
conda env create -f MD-IFP.yml


## Scripts:

 __Trajectories.py__  - functions for building a trajectory object for reading and analysis of standard MD and RAMD trajectories and computation of relative residence times

__IFP_generation.py__  -  functions for generation of IFPs

__Clustering.py__   - functions for analysis of trajectories using IFP data   (is still under developments)

__Process_pdb.py__   - preprocessing PDB files (splitting into ligand and protein files)

__chimera_hydrogen_mol2.py__  - generation of ligand mol2 file 

__IFP_preprocess_Gromacs.py__  - enables wrapping a system back into the original box using trjconv Gromacs tool. Script is designed for a specific file structure - please adjust accordingly. The script helps to transform system back into the box in the most  but not in 100% of cases. For example it does not prevent splitting two proteins in the case of protein-protein complexes 

       
## Application examples (folder Examples):

   1. IFP.py - Generation of the IFP databease for a single MD trajectory of a protein-ligand complex
   2. IFP_contacts_quickView.py  - generation of a plot with average IFPs extracted from a trajectory

## Test Examples as Python Jupyter Notebooks :

### Data employed in test examples 
   can be downloaded from  https://zenodo.org/record/3981155#.XzQEUCgzaUk

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

*This plot illustrates ligand dissociation pathways in a graph representation derived from clustering ligand trajectories in IFP space and plotting them with respect to the ligand COM from the initial bound position.*
   
### IV. __MD-IFP_test.ipynb:__

JN designed for validation of the IFP sctipt on 40 PDB complexes (used in paper J. Chem. Phys. 2020)

![HSP90](/images/IFP_validation.png)

*Validation of PL IFP on 40 PDB structures:*
   - false positives (FP) if no other method (FLIP, PLIP, LPC, MOE) was able to find them 
   - false negative (FN) if all four (three for water bridges and halogen bonds) found the missing interaction.
   - true positives (TP) if at least one method (FLIP, PLIP, LPC, MOE) was able to find it
