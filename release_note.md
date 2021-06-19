##  MD-IFP : release note


__version 1.0__   (Kokh et al. JCP 2020) 
__version 1.1__
   patches:
   1. procedure for detection of aromatic interactions revised (several bugs fixed):
   - CG* aromatic atoms missed 
   - cation-aromatic and aromatic-aromatic were mixed out
   - prot.-ligand and ligand - protein dist. 5.5 A
   - ligand S - aromatic residue interaction removed
   - in adding contacts thet were not counted - (< 4) updated to (<= 4)
   2. added JN for testing IF on 40 PDB complexes  
   - selection of atoms by type (for example "type H") was replaced by by selection by name ("name H*")
__version 1.2__ (currently master branch)
   1. Added IFP_preprocess_Gromacs.py 
   2. 
