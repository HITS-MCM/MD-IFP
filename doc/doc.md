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


# __Important notes regarding the file preparation procedure for IFP computations:__

## __Ligand structure preparation__
  1. Ligand should be protonated. In the case of multiple ligands, one can use automated procedure implemented in the script Process_pdb.py that employes Chimera software (https://www.cgl.ucsf.edu/chimera/). 
  2. Ideally ligand structure should be provided by MOL2 file. However, not all MOL2 formats are accepted by RDKit (Python library that is used to determine ligand atom  properties). The best way to generate mol2 file is to use MOE or Maestro software. Generated mol2 file in same cases can also be corrected by http://www.swissparam.ch/. The main problem with Chimera is that correctly describes bonds in aromatic or cyclic groubs containing nitrogen atoms
  3. If mol2 file is absent or is not accepted by RDKit, pdb file will be used to define properies of ligand atoms. Unfortunately, in this case aromatic fragments will not be identified.
  4. There are several atom properties that are introduced in addition to those identified by RDKit:
      - Hydrophobic properties of fluorine atoms
      -  Neg. ionazable property of the phosphate atom
      -  Acceptor property of the oxigen atoms bound to phosphate atoms
  5. __Important__, the name of the ligand in mol2 file (or the residue name in the pdb file if mol2 is absent) is used to detect ligand in a trajectory or in a complex. If the residue name is different - IFP will not be computed!
