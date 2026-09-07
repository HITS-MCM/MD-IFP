#!/usr/bin/env python
# coding: utf-8

#---------------HOW TO USE--------------------
# within the folder contaning all the respective mutant folders (e.g., BN-BS/mutants/) execute for example:
#
# python IFP_PP-eq.py WT > IFP-EQ_WT.log ; where
#
# "WT" is sys.argv[1] corresponding to the WT folder 
#
###############################################################################################################
# 
# Package for analysis of GROMACS equilibration trajectories of a protein-protein complex
# 
################################################################################################################
# 
# Scripts for the generation of IFPs from equilibration trajectories performed with Gromacs.
# The file includes scripts for:
#  - post-processing GROMACS equilibration trajectories by e.g., solving PBC.
#  - Protein-Protein REsidue pair contacts (PP-REs) computaton
#  - RMSD calculation
#  - water analysis
# Ouputs generated:
#  - COM_distance-"+Trj+".dat; for every equilibration trajectory a file with protein-protein COM-COM for every frame
#  - RMSD-"+Trj+".dat; for every equilibration trajectory a file with RMSD
#  - Water Analysis was removed in this version, in order to increase Performance
#
# For usage, please adjust directory/file structure accordingly where indicated.
#
################################################################################################################

# Scripts were used for analyses reported in the manuscript
# "Computation of the effects of mutations on protein-protein dissociation rates and unbinding mechanisms" 
# by Giulia D'Arrigo, Daria B. Kokh, Ariane Nunes-Alves and Rebecca C. Wade
# submitted to *** , 2024

################################################################################################################

# Packages required:
# numpy 
# pandas
# MDAnalysis 
# code is written on Python 3.x and tested on the version 3.7

################################################################################################################

# 02.06.2021
# Copyright (c) 2024
# Released under the EUPL Licence, v1.2 or any higher version

################################################################################################################

# Author: Daria B. Kokh, Giulia D'Arrigo
# mcmsoft@h-its.org
# Heidelberg Institute of Theoretical Studies (HITS, www.h-its.org)
# Schloss-Wolfsbrunnenweg 35
# 69118 Heidelberg, Germany

################################################################################################################

import glob, os
import sys
import subprocess
import numpy as np
from numpy import linalg as LA
from numpy import asarray
from numpy import save


import pandas as pd
from pandas import ExcelFile 
import MDAnalysis as mda
from MDAnalysis.analysis import contacts,align,rms,distances
from MDAnalysis.analysis.base import AnalysisFromFunction
from MDAnalysis.coordinates.memory import MemoryReader
import MDAnalysis.analysis.hbonds as hb
from scipy.spatial import distance_matrix
import MDAnalysis.analysis.rms
from MDAnalysis import transformations
from MDAnalysis.lib.mdamath import make_whole
import matplotlib.pyplot as plt


at_aromatic = "((resname PHE TRP TYR HIS HIE HID HE2) and (name CZ* CD* CE* CG* CH* NE* ND*))"
at_positive =  "((resname ARG LYS ) and (name NH* NZ)) or ((resname HI2 ) and (name HD HE))"
at_negative = " ((resname ASP GLU) and (name OE* OD*))"
at_Hdon = "(resname TYR SER LYS GLN ARG HIS ASN THR and (not backbone) and (name O* N*))"
at_Hacc = " (resname GLU ASP GLN and (not backbone) and (name O*))"
at_sulfur = "(protein and (name S*))"
at_hydrophob = " (protein  and (name C*  S*) and (not  (name CG and resname ASN ASP))   and (not  (name CD and resname GLU GLN ARG))  and (not  (name CZ and resname TYR ARG))  and (not  (name CE and resname LYS)) and (not  (name CB and resname SER THR))   and (not backbone))"
at_BB_positive = " backbone and name N "
at_BB_negative = " backbone and name O "

mut_name  = sys.argv[1]
trj_list_n  = sys.argv[2]
replica_dir = f"Replica-{trj_list_n}/"

####modify accordingly####
PRJ_DIR =  "/hits/fast/mcm/szablods/MD-IFP/RAMD-700-WT-BN-BS"      #insert the directory path containing all the mutants, "mut_name" will be specified in the bash file for running this script
pdb = "/hits/fast/mcm/szablods/MD-IFP/BN-BS/mutants/WT/ref_with_segids.pdb"                 #insert the directory path containing the ref.pdb structure (prepared with tleap),
top = "/hits/fast/mcm/szablods/MD-IFP/RAMD-700-WT-BN-BS/Gromacs-large/gromacs.top"             #insert the directory path containing the gromacs topology

#create a folder for EQ outputs
OUTPUT_DIR = "/hits/fast/mcm/szablods/MD-IFP/"
os.makedirs(OUTPUT_DIR, exist_ok=True)

#resi_lig = [1,109]                   #insert the residue index range for protein 1 (the one where the force is applied, named as lig)
#resi_target = [110,198]              #insert the residue index range for protein 2 (named as target)

sel_lig = "segid PROA and (not name H*)"               # your ligand
sel_target = "segid PROB and (not name H*)"  # receptor chains
selection = "(segid PROA or segid PROB) and not name H*"  # for whole system reduced
selection1 =  " not name H*"

ramd_dir = "/RAMD-700/"          #insert the directory path containing the GROMACS trajectories of interest
ramd_trj = "traj_comp.xtc"
ramd_tpr = "gromacs_ramd.tpr"
ramd_trj_fixed1 = "tr1.xtc"
ramd_trj_fixed2 = "tr2.xtc"
ramd_trj_fixed = "tr.xtc"
ramd_trj_fixed0 = "tr0.xtc"
Trj_list = []
for nn in range(1,20): Trj_list.append("Replica"+str(nn))   ##number of GROMACS equilibration replicas to analyse

Trj_list = {"TRJ-1-1", "TRJ-1-2", "TRJ-1-3", "TRJ-1-4", "TRJ-1-5", "TRJ-1-6", "TRJ-1-7", "TRJ-1-8", "TRJ-1-9", "TRJ-1-10", "TRJ-2-1", "TRJ-2-2", "TRJ-2-3", "TRJ-2-4", "TRJ-2-5", "TRJ-2-6", "TRJ-2-7", "TRJ-2-8", "TRJ-2-9", "TRJ-2-10", "TRJ-3-1", "TRJ-3-2", "TRJ-3-3", "TRJ-3-4", "TRJ-3-5", "TRJ-3-6", "TRJ-3-7", "TRJ-3-8", "TRJ-3-9", "TRJ-3-10"}

for Trj in Trj_list :
  print("============================================================")
  print(Trj,ramd_dir+Trj+"/"+ramd_trj)
  trj = ramd_dir+Trj+"/"+ramd_trj
  if not os.path.exists(PRJ_DIR+trj):
      print("trajectory was not found ", PRJ_DIR+trj)
      continue
  else:
      print(">>>>>>>>>>>>>>>>>>>>>>",trj)
  if not os.path.exists(PRJ_DIR+ramd_dir+Trj+"/"+ramd_trj_fixed):
    with open(PRJ_DIR+ramd_dir+Trj+"/preprocess.sh", 'w') as f_bash:
       f_bash.write("#!/bin/bash\n")
       f_bash.write("cd "+PRJ_DIR+ramd_dir+Trj+" ;  gmx trjconv -f "+ramd_trj+"  -s "+ramd_tpr +"  -pbc res -o "+ramd_trj_fixed1 +"  <<< 0\n")
       f_bash.write("gmx trjconv -f "+ramd_trj_fixed1+"  -s "+ramd_tpr +"  -pbc nojump -o "+ramd_trj_fixed2 +" <<< 0\n")
       f_bash.write("gmx trjconv -f "+ramd_trj_fixed2+"  -s "+ramd_tpr +"  -pbc whole -o "+ramd_trj_fixed +" <<< 0\n")
     #  f_bash.write("gmx trjconv -f "+ramd_trj+"  -s "+ramd_tpr +"  -pbc nojump -o "+ramd_trj_fixed3 +" <<< 0\n")
       f_bash.write("gmx trjconv -f "+ramd_trj_fixed2+"  -s "+ramd_tpr +"  -pbc atom -o "+ramd_trj_fixed0 +" <<< 0\n")
       f_bash.write("rm "+ramd_trj_fixed1+" "+ramd_trj_fixed2)
    os.chmod(PRJ_DIR+ramd_dir+Trj+"/preprocess.sh" , 0o755)
    os.system(PRJ_DIR+ramd_dir+Trj+"/preprocess.sh")
  trj = ramd_dir+Trj+"/"+ramd_trj_fixed
  trj_wat = ramd_dir+Trj+"/"+ramd_trj_fixed0

  ref = mda.Universe(pdb)
  u = mda.Universe(pdb,PRJ_DIR+trj,guess_bonds=False)


###### COM-COM distance and PP-REs calculation


  system_reduced = u.select_atoms(selection)
  start = 0
  step = 1
  step_beg = 2
  stop = len(u.trajectory)

  u_length = len(u.trajectory)
  u_size = int(os.path.getsize(PRJ_DIR+trj)/(1024.*1024.))
  print("total number of frames= %s; file size %s M" %(u_length,u_size))
  u_mem = mda.Merge(system_reduced).load_new(AnalysisFromFunction(lambda ag: ag.positions.copy(), system_reduced).run(start=start,stop=stop,step=step).results,format=MemoryReader)
  u = u_mem
  dist_file = open(os.path.join(OUTPUT_DIR, "COM_distance-"+Trj+".dat"), 'w')
  distance_COM_list = []
  matr_dist = []
  times_saved = []

  print(".....preparation...")
  u_CA_lig_IP = []
  u_CA_lig_IN = []
  u_CA_lig_HD = []
  u_CA_lig_HA = []
  u_CA_lig_AR = []
  u_CA_lig_HY = []
  u_CA_lig_BB = []
  u_CA_lig_BB_P = []
  u_CA_lig_BB_N = []
  u_CA_lig_C = []
  u_CA_target_IP = []
  u_CA_target_IN = []
  u_CA_target_HD = []
  u_CA_target_HA = []
  u_CA_target_AR = []
  u_CA_target_HY = []
  u_CA_target_BB = []
  u_CA_target_BB_P = []
  u_CA_target_BB_N = []
  u_CA_target_C = []

  for res in u.select_atoms("segid PROC").residues:
    u_CA_lig_IP.append(u.select_atoms(f"resid {res.resid} and {at_positive}", updating=False))
    u_CA_lig_IN.append(u.select_atoms(f"resid {res.resid} and {at_negative}", updating=False))
    u_CA_lig_HD.append(u.select_atoms(f"resid {res.resid} and {at_Hdon}", updating=False))
    u_CA_lig_HA.append(u.select_atoms(f"resid {res.resid} and {at_Hacc}", updating=False))
    u_CA_lig_AR.append(u.select_atoms(f"resid {res.resid} and {at_aromatic}", updating=False))
    u_CA_lig_HY.append(u.select_atoms(f"resid {res.resid} and {at_hydrophob}", updating=False))
    u_CA_lig_BB.append(u.select_atoms(f"resid {res.resid} and backbone", updating=False))
    u_CA_lig_BB_N.append(u.select_atoms(f"resid {res.resid} and backbone and name O*", updating=False))
    u_CA_lig_BB_P.append(u.select_atoms(f"resid {res.resid} and backbone and name N*", updating=False))
    u_CA_lig_C.append(u.select_atoms(f"resid {res.resid} and name C* and not backbone", updating=False))

  for res in u.select_atoms("(segid PROA or segid PROB)").residues:
      u_CA_target_IP.append(u.select_atoms(f"resid {res.resid} and {at_positive}", updating=False))
      u_CA_target_IN.append(u.select_atoms(f"resid {res.resid} and {at_negative}", updating=False))
      u_CA_target_HD.append(u.select_atoms(f"resid {res.resid} and {at_Hdon}", updating=False))
      u_CA_target_HA.append(u.select_atoms(f"resid {res.resid} and {at_Hacc}", updating=False))
      u_CA_target_AR.append(u.select_atoms(f"resid {res.resid} and {at_aromatic}", updating=False))
      u_CA_target_HY.append(u.select_atoms(f"resid {res.resid} and {at_hydrophob}", updating=False))
      u_CA_target_BB.append(u.select_atoms(f"resid {res.resid} and backbone", updating=False))
      u_CA_target_BB_N.append(u.select_atoms(f"resid {res.resid} and backbone and name O*", updating=False))
      u_CA_target_BB_P.append(u.select_atoms(f"resid {res.resid} and backbone and name N*", updating=False))
      u_CA_target_C.append(u.select_atoms(f"resid {res.resid} and name C* and not backbone", updating=False))


  for i in range(0,len(u.trajectory)):
    u.trajectory[i]
    if (i < len(u.trajectory)/2) and i%step_beg > 0: continue
    if (i < len(u.trajectory)/2): interval = step_beg
    else: interval = step
    u.dimensions = u.dimensions
    u_CA_lig = u.select_atoms(sel_lig )
    u_CA_target = u.select_atoms(sel_target )

    distance_COM = LA.norm (u_CA_lig.center_of_mass()- u_CA_target.center_of_mass())
    print("DISTANCE = ", i,distance_COM)
    dist_file.write("%i %6.1f %i \n" %(i,distance_COM,interval))
    distance_COM_list.append(distance_COM)
    matr = []
    for l in range(0, len(u_CA_lig_IP)):
        lt = []
        for t in range(0, len(u_CA_target_IP)):
             # possible polar bond
             if len(u_CA_target_C[t]) == 0 or len(u_CA_lig_C[l]) == 0:
                  d = int(10* LA.norm (u_CA_lig_BB[l].center_of_mass()- u_CA_target_BB[t].center_of_mass()))
             else:
                  d = int(10* LA.norm (u_CA_lig_C[l].center_of_mass()- u_CA_target_C[t].center_of_mass()))
             tt = []
             tt.append(str(d)+"--")
             if d < 150:
               if len(u_CA_lig_IP[l]) > 0 and len(u_CA_target_IN[t]) > 0:
                  d1 = int(10* LA.norm (u_CA_lig_IP[l].center_of_mass()- u_CA_target_IN[t].center_of_mass()))
                  tt.append("IP-IN"+str(d)+"-"+str(d1))
                  d = min(d1,d)
               if len(u_CA_lig_IN[l]) > 0 and len(u_CA_target_IP[t]) > 0:
                  d1 = int(10* LA.norm(u_CA_lig_IN[l].center_of_mass()- u_CA_target_IP[t].center_of_mass()))
                  tt.append("IN-IP"+str(d)+"-"+str(d1))
                  d = min(d1,d)
               if len(u_CA_lig_HD[l]) > 0 and len(u_CA_target_HA[t]) > 0:
                  d1 = int(10* LA.norm(u_CA_lig_HD[l].center_of_mass()- u_CA_target_HA[t].center_of_mass()))
                  tt.append("HD-HA"+str(d)+"-"+str(d1))
                  d = min(d1,d)
               if len(u_CA_lig_HA[l]) > 0 and len(u_CA_target_HD[t]) > 0:
                  d1 = int(10* LA.norm(u_CA_lig_HA[l].center_of_mass()- u_CA_target_HD[t].center_of_mass()))
                  tt.append("HA-HD"+str(d)+"-"+str(d1))
                  d = min(d1,d)
               if len(u_CA_lig_AR[l]) > 0 and len(u_CA_target_AR[t]) > 0:
                  d1 = int(10* LA.norm(u_CA_lig_AR[l].center_of_mass()- u_CA_target_AR[t].center_of_mass()))
                  tt.append("AR-AR"+str(d)+"-"+str(d1))
                  d = min(d1,d)
               if len(u_CA_lig_HY[l]) > 0 and len(u_CA_target_HY[t]) > 0:
                  d1 = int(10* LA.norm(u_CA_lig_HY[l].center_of_mass()- u_CA_target_HY[t].center_of_mass()))
                  tt.append("HY-HY"+str(d)+"-"+str(d1))
                  d = min(d1,d)
               if len(u_CA_lig_BB_P[l]) > 0 and len(u_CA_target_HA[t]) > 0:
                  d1 = int(10* LA.norm(u_CA_lig_BB_P[l].center_of_mass()- u_CA_target_HA[t].center_of_mass()))
                  tt.append("BB-HA"+str(d)+"-"+str(d1))
                  d = min(d1,d)
               if len(u_CA_lig_BB_N[l]) > 0 and len(u_CA_target_HD[t]) > 0:
                  d1 = int(10* LA.norm(u_CA_lig_BB_N[l].center_of_mass()- u_CA_target_HD[t].center_of_mass()))
                  tt.append("BB-HD"+str(d)+"-"+str(d1))
                  d = min(d1,d)
               if len(u_CA_lig_HD[l]) > 0 and len(u_CA_target_BB_N[t]) > 0:
                  d1 = int(10* LA.norm(u_CA_lig_HD[l].center_of_mass()- u_CA_target_BB_N[t].center_of_mass()))
                  tt.append("HD-BB"+str(d)+"-"+str(d1))
                  d = min(d1,d)
               if len(u_CA_lig_HA[l]) > 0 and len(u_CA_target_BB_P[t]) > 0:
                  d1 = int(10* LA.norm(u_CA_lig_HA[l].center_of_mass()- u_CA_target_BB_P[t].center_of_mass()))
                  tt.append("HA-BB"+str(d)+"-"+str(d1))
                  d = min(d1,d)
             if d < 35: print("FOUND:", l, t, tt, d)             
             lt.append(d)
          #   if d < 50:
          #       print(l,t,"--",d,"--",len(u_CA_lig_NO),len(u_CA_target_NO), len(u_CA_lig_BB), len(u_CA_target_BB) , len(u_CA_lig_C), len(u_CA_target_C))
        matr.append(lt)
    matr_dist.append(matr)
    times_saved.append(i)
  np.save(os.path.join(OUTPUT_DIR, "matrix_distance-"+Trj+".npy"), np.asarray(matr_dist).astype(dtype='int16'))
  dist_file.close()
  print(np.asarray(matr_dist).astype(dtype ='int16').shape, len(distance_COM_list))
###### RMSD

  R = MDAnalysis.analysis.rms.RMSD(u, ref,
           select=" backbone and "+sel_target,             # superimpose on whole backbone of the whole protein
           groupselections=["(not name H*) and "+sel_lig,   # CORE
                            " (not name H*) and "+sel_target])                                    # NMP
  R.run()
  rmsd = R.rmsd.T
  R1 = MDAnalysis.analysis.rms.RMSD(u, ref,
           select=" backbone and "+sel_lig,             # superimpose on whole backbone of the whole protein
           groupselections=["(not name H*) and "+sel_lig,   # CORE
                            " (not name H*) and "+sel_target])                                    # NMP
  R1.run()
  rmsd1 = R1.rmsd.T

  time = rmsd[1]
  fig = plt.figure(figsize=(4,4),dpi=150)
  ax = fig.add_subplot(111)
#ax.plot(time, rmsd[2], 'b-',  label="all")
  lig_label = "PROC"  # your ligand segment
  target_label = "PROA+PROB"  # your receptor segments
  ax.plot(time, rmsd[3], 'k-',  label=f"segid-{target_label}")
  ax.plot(time, rmsd[4], 'k--', label=f"segid-{target_label}")
  ax.plot(time, rmsd1[3], 'b-',  label=f"segid-{lig_label}")
  ax.plot(time, rmsd1[4], 'b--', label=f"segid-{lig_label}")
  ax.plot(time[np.array(times_saved)],distance_COM_list, 'r--', label="COM-COM")
  ax.legend(loc="best")
  ax.set_ylim(0,70)
  ax.set_xlabel("time (ps)")
  ax.set_ylabel("RMSD ($\AA$)")
  fig.savefig(os.path.join(OUTPUT_DIR, Trj+".jpg"))
  
  rmsd_file = open(os.path.join(OUTPUT_DIR, "RMSD-"+Trj+".dat"), 'w')
  for r0, r1,r2, r3, r4, r5, r6 in zip(time, rmsd1[2], rmsd1[3], rmsd1[4], rmsd[2] , rmsd[3],  rmsd[4]):
     rmsd_file.write("%6.1f %6.1f %6.1f %6.1f %6.1f %6.1f  %6.1f\n" %(r0,r1,r2,r3,r4,r5,r6))
  rmsd_file.close()
