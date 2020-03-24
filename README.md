# MD-trajectory_analysis
A workflow for generation and analysis Protein-Ligand Interaction Fingerprints from MD tajectories


    Author: Daria Kokh
    Test examples were revised by: Fabian Ormersbach 


    v.1.0
    Copyright (c) 2019
    Released under the GNU Public Licence, v2 or any higher version
    Daria.Kokh@h-its.org
    Heidelberg Institute of Theoretical Studies (HITS, www.h-its.org)
    Schloss-Wolfsbrunnenweg 35
    69118 Heidelberg, Germany
    

This open source software code was developed in part in the Human Brain Project, funded from the European Union’s Horizon 2020 Framework Programme for Research and Innovation under Specific Grant Agreements  No. 785907 (Human Brain Project  SGA2).

    Packages required:
    Python 3.x
    numpy;    pandas;  matplotlib;  MDAnalysis;  RDkit;   scipy
    code is written on Python 3.x and tested on Python 3.7
    to configure envirement in anaconda use
    conda env create -f IFP_trajectory.yml



1. Scripts:

    Clustering.py   - containes functions for analysis of trajectories using IFP data   
    IFP_generation.py  - contains functions for generation of IFPs
    Membrane_analysis.py
    Trajectories.py  - contains functions for building a trajectory object for reading and analysis standard MD and RAMD trajectories and computation of relative residence times

2. Data:

        2YKI
        4MQT
        6EI5
        
3. Example jupyter notebooks :

         (i)jupyter notebook IFP_generation_examples_PDB.ipynb: an examples of Protein-Ligand IFP computations
        Computing interaction fingerprints (IFP) for
        --- a single structure prepared for MD simulations (HSP90; PDB ID 6EI5, dcd format)
        --- a trajectory (for selected frames; dcd format)
        --- a PDB structure

        (ii)jupyter notebook IFP_generation_example_TRAJ.ipynb: example of generation and analysis of Ligand-Protein IFPs for RAMD simulations (dcd trajectories)
        Computing interaction fingerprints (IFP) for
        --- system equilibration trajectory
        --- ligand dissociation trajectories
