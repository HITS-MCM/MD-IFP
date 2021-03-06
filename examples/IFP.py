#!/usr/bin/env python
# coding: utf-8

#  Script for generation of protein-ligand interaction fingerprints, IFPs, from MD trajectories
# 
# 
#############################
### v 1.0
#
#    Copyright (c) 2020
#    Released under the GNU Public Licence, v2 or any higher version
#    
### Author: Daria Kokh
#    Daria.Kokh@h-its.org
#    Heidelberg Institute of Theoretical Studies (HITS, www.h-its.org)
#    Schloss-Wolfsbrunnenweg 35
#    69118 Heidelberg, Germany
################################# 
# 
# ### Input data required:
#     trajectory file (possible formats: nc,dcd, xtc, trr )
#     pdb file of the system (for example, generated from the first frame)
#     ligand pdb file
#     ligand mol2 file (not all mol2 files can be read, files generated by MOEor Maestro are fine)
#     

import glob, os, sys
from Scripts.IFP_generation import *
from Scripts.Trajectories import *


####
# Location of your data:
#  your_path/ligand_name/top/input.pdb - reference PDB structure
#  your_path/ligand_name/build/moe.mol2  - ligand mol2 file
#  your_path/ligand_name/build/ligand_name.mol2  - ligand pdb file
#  your_path/ligand_name/md/out/*.nc    - trajectories to be analyzed
###



ligand_name = name
DIR_all = "your_path"+name+"/"
trj = "md/out/"
ref_pdb = "top/input.pdb"
ligand_pdb = "build/"+name+".pdb"
ligand_mol2 = "build/moe.mol2"

#os.system("grep "+ligand_name+" "+ref_pdb+ " > build/"+name+"_resp.pdb")

step =1
start=0
end = -1
print("TRJ "+name+">>>>"+DIR_all+trj)
sys.stdout.flush()
tr = trajectories(DIR_all,namd_tmpl= trj, ramd_tmpl= DIR_all,ligand_pdb=ligand_pdb,pdb = ref_pdb,\
                          ligand_mol2=ligand_mol2,namd_traj_tmpl = "*nc",ramd_traj_tmpl = "xx")
tr.sub_system = " protein or (resname SOL HOH WAT G G3 G5 U5 C C3 MN) "
print("TRJ "+name+">> IFP")
sys.stdout.flush()
tr.analysis_all_namd(WB_analysis = False, Lipids = [],auxi_selection = [],step_analysis=step, start_analysis=start)

sys.stdout.flush()
IFP_table = tr.namd.IFP_save(DIR_all+ligand_name+"-IFP_3000.pkl")

