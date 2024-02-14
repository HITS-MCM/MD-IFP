# MD-IFP: Interaction Fingerprint analysis of Molecular Dynamics trajectories
<br><br>

## Description
The MD-IFP is a python workflow for the generation and analysis of protein-ligand (PL) and protein-protein (PP) interaction fingerprints from Molecular Dynamics trajectories. If used for the analysis of RAMD (Random Accelaration Molecular Dynamics) trajectories, it can help to investigate dissociation mechanisms by characterizing transition states as well as the determinants and hot-spots for dissociation. As such, the combined use of τRAMD and MD-IFP may assist the early stages of drug discovery campaigns for the design of new molecules or ligand optimization.
<br><br>

## Cite us
If you use or adapt MD-IFP for your own research projects please cite us.
*D. B. Kokh, B. Doser, S. Richter, F. Ormersbach, X. Cheng , R.C. Wade  "A Workflow for Exploring Ligand Dissociation from a Macromolecule: Efficient Random Acceleration Molecular Dynamics Simulation and Interaction Fingerprints Analysis of Ligand Trajectories" J. Chem Phys.(2020) 158  125102  doi: 10.1063/5.0019088.*
<br><br>  

## Publications of the method application: 
1. IFP analysis of dissociation trajectories for 3 compounds of HSP90  reported in the paper 
  
   D. B. Kokh, B. Doser, S. Richter, F. Ormersbach, X. Cheng , R.C. Wade  "A Workflow for Exploring Ligand Dissociation from a Macromolecule: Efficient Random Acceleration Molecular Dynamics Simulation and Interaction Fingerprints Analysis of Ligand Trajectories" J. Chem Phys.(2020) 153  125102  doi: 10.1063/5.0019088; https://arxiv.org/abs/2006.11066
   
   Results are implemented in [Dissociation mechanisms of HSP90-small molecules](./Protein-Ligand/Examples-JN/IFP_analysis_RAMD_HSP90.ipynb)
  
2.  Small compound unbinding from T4 lysozyme mutants

    A Nunes-Alves, DB Kokh, RC Wade  "Ligand unbinding mechanisms and kinetics for T4 lysozyme mutants from τRAMD simulations", Current Research in Structural Biology 3, 106-111
    https://doi.org/10.1016/j.crstbi.2021.04.001

3. Application to two GPCR targets (embedded in a membrane):

   D. B. Kokh, R.C. Wade "G-Protein Coupled Receptor-Ligand Dissociation Rates and Mechanisms from τRAMD Simulations", doi: https://doi.org/10.1101/2021.06.20.449151

   Associated scripts and data can be downloaded at https://zenodo.org/record/5001884#.YM-rRmgzYuU
   
   Results are implemented in [Dissociation mechanisms from M2 and β2 adrenergic receptors](./Protein-Ligand/Examples-JN/IFP_analysis_RAMD_GPCR.ipynb)

4. Application to Protein-Protein complexes:
   
   Results are implemented in [Dissociation mechanisms of protein-protein complexes](./Protein-Protein/Examplex-JN/IFP_analysis_RAMD_PP.ipynb)
<br><br>

## Packages requirements:
__Python 3.x__ 

__Python Libraries:__ 
   1. numpy; pandas; matplotlib; seaborn; sklearn; scipy; 
   2. __RDkit__ - only for the MD-IFP analysis of protein-ligand complexes
   3. __ngview__ - used for visualization. Installation of ngview can be tricky, the following way may work: after installation of the Python envirenment, execute:
```
conda install -c conda-forge nglview=2.7.1
```
and then

```
jupyter-nbextension enable nglview --py --sys-prefix
```
If you don't need visualization, you can skip this, but JN must be edited accordingly.

   4. __MDAnalysis__ Version 1.1.1 - (Important: an old module for H-bond analysis is currently used, it will be removed in version 2.0 ). Best is to use the __MD-IFP.yml__ file to generate a python environment in anaconda (as shown below).

__Chimera__ - only for the scripts used for preprocessing pdb files (structure protonation and generation of the ligand mol2 file; not required if protonation and mol2 file are already prepared by a user)
    
__Codes were written on Python 3.x and tested on Python 3.8__

__To configure environment in anaconda use:__

```
conda env create -f MD-IFP.yml
```
<br><br>

## Acknowledgments

This open source software code was developed in part in the __Human Brain Project__, funded from the European Union’s Horizon 2020 Framework Programme for Research and Innovation under Specific Grant Agreements  No. 785907 (Human Brain Project  SGA2).
