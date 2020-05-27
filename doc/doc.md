# MD-IFP: MD trajectory analysis using protein-ligand Interaction Fingerprints
## A Python workflow for generation and analysis Protein-Ligand Interaction Fingerprints from Molecular Dynamics tajectories
## v.1.0

__Authors:__

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

__To configure envirement in anaconda use:__
conda env create -f IFP_trajectory.yml


##__Important Notes:__
