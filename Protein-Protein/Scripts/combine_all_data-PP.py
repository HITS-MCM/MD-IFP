#!/usr/bin/env python
# coding: utf-8


################################################################################################################
# 
# Script for the generation of .pkl files of equilibration and RAMD trajectories, generated with IFP_PP-EQ.py and IFP_PP.py respectively, for each analysed mutant.
# Ouputs generated:
#  - "Mutant"-EQ.pkl and "Mutant"-EQ_water.pkl; for every mutant, equilibration pkl files
#  - "Mutant"-RAMD.pkl and "Mutant"-RAMD_water.pkl; for every mutant, RAMD pkl files
#
# For usage, please adjust the fields where indicated.
# Pay attention to adjust the residue numbering in "resi_update" when the "read_repl_trj" function is called for both equilibration and RAMD,
# in this way the residue numbering will correspond to the original crystal structure (to ease understanding of IFP results).
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

import pandas as pd
from pandas import ExcelFile 



def read_repl_trj(dir_all, threshold=150,resi_update=[0,0],RAMD=True):
    """
    Parameters:
    Results:
    """
    print(dir_all)
    COM = []
    RMSD = []
    time = []
    interval = []
    timeW = []
    W = []
    Wt = []
    trW = []
    for r in range(0,10):
        for t in range(0,20):
            if RAMD:
                suffix = "TRJ-"+str(r+1)+"-"+str(t+1)
            else:
                suffix ="Replica"+str(r+1)
            if (not os.path.exists(dir_all+"matrix_distance-"+suffix+".npy")) or (not os.path.exists(dir_all+"Water-"+suffix+".dat")): 
                continue
            if os.stat(dir_all+"COM_distance-"+suffix+".dat").st_size == 0 or os.stat(dir_all+"Water-"+suffix+".dat").st_size == 0:
                continue
            with open(dir_all+"COM_distance-"+suffix+".dat", 'r') as f:
                lines = f.readlines()
                time_local = []
                for l in lines: 
                    COM.append(float(l.split()[1]))
                    time.append(int(l.split()[0]))
                    time_local.append(int(l.split()[0]))
                    try:
                        interval.append(int(l.split()[2]))
                    except:
                        interval.append(1)
                if len(time_local)== 0:
                    print (" File Empty- ",r+1,t+1)
                    continue
            with open(dir_all+"Water-"+suffix+".dat", 'r') as f:
                lines = f.readlines()
                for l in lines:  
                    Wt.append(int(l.split()[1]))
                    W.append(int(l.split()[2]))
                    timeW.append(int(l.split()[0]))
                    trW.append(str(r+1)+"-"+str((r+1)*100+t+1))
            s = np.load(dir_all+"matrix_distance-"+suffix+".npy")
            shape1 = s.shape[1]
            shape2 = s.shape[2]
            s = np.reshape(s, (s.shape[0],s.shape[1]*s.shape[2]))
            if len(RMSD) == 0: 
                print(r,t,"START")
                matrix = s
                tr = np.full((s.shape[0]), str(r+1)+"-"+str((r+1)*100+t+1))
            else: 
                matrix = np.append(matrix,s,axis=0)
                tr = np.append(tr,np.full((s.shape[0]), str(r+1)+"-"+str((r+1)*100+t+1)))
            with open(dir_all+"RMSD-"+suffix+".dat", 'r') as f:
                rmsd = []
                for l in f.readlines(): 
                    rmsd.append(l.split()[2]+"-"+l.split()[3]+"-"+l.split()[4]+"-"+l.split()[5])
                rmsd = np.array(rmsd)[np.array(time_local)]
                if  len(RMSD) == 0: print(r,t,"START")
                RMSD = np.append(RMSD,rmsd)
            if not RAMD : break
    if len(RMSD) == 0: 
        print(r,"Files were not found in ",dir_all)
        sys.exit()
    Mtx = np.clip(matrix, 0, threshold, matrix)
    cols_indx = []
    for i in range(0,shape1):
        for j in range(0,shape2):
            cols_indx.append((i+1+resi_update[0],j+1+resi_update[1]))
    df = pd.DataFrame(data=Mtx, columns=cols_indx)
    print(df.shape,len(COM),len(RMSD))
    df['time'] = np.array(time)
    df['Repl'] = np.array(tr)
    df['Traj'] = df['Repl'].str.split('-', 1).str[1]
    df['Repl'] = df['Repl'].str.split('-', 1).str[0]
    df['COM'] = COM
    df['Interval'] = interval
    df['RMSD'] = RMSD
    df['RMSD11'] = df['RMSD'].str.split('-', 1).str[0].astype(float)
    df['RMSD12'] = df['RMSD'].str.split('-', 2).str[1].astype(float)
    df['RMSD21'] = df['RMSD'].str.split('-', 3).str[2].astype(float)
    df['RMSD22'] = df['RMSD'].str.split('-', 4).str[3].astype(float)

    df['min'] = np.min(df[cols_indx].values, axis=1)
    df = df.drop(columns=['RMSD'])
    
    dfW = pd.DataFrame(data=np.array(timeW), columns=["time"])  
    dfW["WaterTot"] = np.array(Wt)
    dfW["Water"] = np.array(W)
    dfW['Repl'] = np.array(trW)
    dfW['Traj'] = dfW['Repl'].str.split('-', 1).str[1]
    dfW['Repl'] = dfW['Repl'].str.split('-', 1).str[0]
    
    return(df,dfW)

dir_all = "/hits/basement/mcm/kokhda/Prot-Prot/BN-BS/"  #insert the directory path containing all the mutant folders

mu_list =  ["WT1", "D35A1","D39A1","H102A1","R59A1", "D54A1", "K27D39A2","K27A2",\
        "E60A1","E80A1","H102Y29A1","H102D39A1","E73A1","K27D39A1","K27A1","K27Y29A1", "K27D35A1","E73A1","H102T42A","Y29A2","T42A1","R87A1","N58A1","K27E76A1"]  #insert the list of mutant to analyse
threshold=150


###### pkl generation of equilibration trajectories
###### please modify residue numbering in "resi_update" for both equilibration and RAMD. 
###### For example resi_update=[1,0] means that in the crystal structure protein1 will start with 2 and protein2 will start with 1
###### and so will be in the generated pkl files.

for mu in mu_list:
   df,dfW= read_repl_trj(dir_all+mu+"/", threshold=threshold,resi_update=[1,0],RAMD=False)   
   df['Pr'] = np.full((df.shape[0]),mu)
   df.to_pickle(dir_all+mu+"-EQ.pkl")
   dfW.to_pickle(dir_all+mu+"-EQ_water.pkl")

print("============================================ RAMD ===========================")
for mu in mu_list:
   df,dfW= read_repl_trj(dir_all+mu+"/", threshold=threshold,resi_update=[1,0],RAMD=True)
   df['Pr'] = np.full((df.shape[0]),mu)
   df_new = pd.DataFrame()
   dfW_new = pd.DataFrame()
   cols = []
   [cols.append(c) for c in df.columns.values if str(c[0]).isdigit()]


   for tr in df.Traj.unique():
        tau = []
        df_r = df[df.Traj == tr]
        dfW_r = dfW[dfW.Traj == tr]
        # skip trajectories without dissociation
        try:
           time_motion = min(0,df_r[df_r.COM > 24.5].time.values[0]-100)
           df_new = df_new.append(df_r[df_r.time > time_motion],sort=True)           
           dfW_new = dfW_new.append(dfW_r,sort=True)
        except: print("No dissociation: ",tr,df_r.COM.values)
   cols_2delete = df_new[cols].columns[(df_new[cols].mean(axis=0) == threshold).values]
   df_new = df_new.drop(columns=cols_2delete )
   
   df_new.to_pickle(dir_all+mu+"-RAMD.pkl")
   dfW.to_pickle(dir_all+mu+"-RAMD_water.pkl")

