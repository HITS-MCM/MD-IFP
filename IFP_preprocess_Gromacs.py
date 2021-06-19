################################################################################################################
# Preprocess of dissociation trajectories generated using Gromacs 
# enables wrapping system back into the box
# scripts are designed for a specific file structure, for usage, please adjust accordingly!
################################################################################################################
# "G-Protein Coupled Receptor-Ligand Dissociation Rates and Mechanisms from tauRAMD Simulations" 
# by Daria B. Kokh, Rebecca C. Wade
# submitte to JCTC , June 2021

################################################################################################################
#Packages require:
#numpy ( 1.18.1)
#pandas (1.0.2)
#RdKit
#code is written on Python 3.x and tested on the version 3.7

################################################################################################################
#02.06.2021
#Copyright (c) 2020
#Released under the EUPL Licence, v1.2 or any higher version

################################################################################################################

#Author: Daria Kokh
#Daria.Kokh@h-its.org or mcmsoft@h-its.org
#Heidelberg Institute of Theoretical Studies (HITS, www.h-its.org)
#Schloss-Wolfsbrunnenweg 35
#69118 Heidelberg, Germany

################################################################################################################
import glob, os
import sys
import subprocess
import numpy as np

DIR_all = "./"
DIR_ramd = "RAMD-2020.2/TRJ1-*"   # a set of dirrectories with RAMD trajectories  traj_comp.xtc

frame2save = 2001    #  the last 2001 frames will be analyzed and saved as  traj_comp_whole.xtc

ramd_trj = "traj_comp.xtc"
ramd_trj_fixed = "traj_comp_whole.xtc"
ramd_tpr = "gromacs_ramd.tpr"
ramd_trj_fixed1 = "traj_comp_whole1.xtc"
ramd_trj_fixed2 = "traj_comp_whole2.xtc"


for dir_trj in glob.glob(DIR_all+DIR_ramd):
    if os.path.isdir(dir_trj):
        start_time = -2
        print("----------------------",dir_trj, glob.glob("*out"))
        os.chdir(dir_trj)
        for out_file in glob.glob("*out"):
            os.system(" grep \" GROMACS will be stopped\"  " +out_file+" > tmp")
        try:
            with open(dir_trj+"/tmp", 'r') as f:
                r = f.readline()
                if len(r) > 0:
                    frames = int(r[r.find("after")+6:r.find("steps")-1])/500 # saved frames
                    start_time = int((frames - frame2save))   # start from in ps
                    print(dir_trj,">>>>>>>>>>>>>>>>>>>>>>>>>>>>>",frames,start_time)
                else:
                    start_time = -1
        except:
            print(dir_trj,">>>>>>>>>>>>>>>>>>>>>>>>>>>>> NO OUT file found")
            start_time = -1
        if start_time == -1:
            print(">>>>>>>>>>>>>>>>>>>>>>>>>>>> Skip ",dir_trj)
            continue
        elif start_time < 0:
            start_time = 0
    with open(dir_trj+"/preprocess.sh", 'w') as f_bash:
       f_bash.write("#!/bin/bash\n")
       f_bash.write("cd "+dir_trj+" ;  gmx trjconv -f "+ramd_trj+"  -s "+ramd_tpr +"  -pbc nojump -o "+ramd_trj_fixed1 +" -b "+str(start_time)+" <<< 0\n")
       f_bash.write("cd "+dir_trj+" ;  gmx trjconv -f "+ramd_trj_fixed1+"  -s "+ramd_tpr +"  -pbc atom -o "+ramd_trj_fixed2 +" <<< 0\n")
       f_bash.write("cd "+dir_trj+" ;  gmx trjconv -f "+ramd_trj_fixed2+"  -s "+ramd_tpr +"  -pbc mol -o "+ramd_trj_fixed +" <<< 0\n")
    os.chmod(dir_trj+"/preprocess.sh" , 0o755)
    os.system(dir_trj+"/preprocess.sh")

