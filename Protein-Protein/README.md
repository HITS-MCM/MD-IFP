# A Python Workflow for the Generation and Analysis of Protein-Protein Interaction Fingerprints from Molecular Dynamics trajectories
## This is a collection of examples showing how to carry out particular tasks using MD-IFP scripts
### v.1.1
### 07.02.2024

## Authors and Contributors:

* Daria Kokh
* Giulia D'Arrigo

Heidelberg Institute of Theoretical Studies (HITS, www.h-its.org)

Schloss-Wolfsbrunnenweg 35

69118 Heidelberg, Germany


## Scripts:

__JN_scripts.py__  - Functions needed to run the Jupyter Notebook

__IFP_PP-EQ.py__  - Python script to post-process equilibration trajectories and compute IFP

__IFP_PP.py__   - Python script to post-process RAMD trajectories and compute IFP

__combine_all_data-PP.py__   - Python script to general equilibration and RAMD .pkl files needed to post-analyse trajecotries with the Jupyter Notebook


## Test Examples as Python Jupyter Notebooks :

### I. __IFP_analysis_RAMD_PP.ipynb:__

JN for post-processing RAMD trajectories and explore dissociation mechanisms of Bn-Bs and BT/BCT-BPTI mutants. This JN uses pre-generated IFP databases.
