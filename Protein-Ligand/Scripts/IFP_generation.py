#!/usr/bin/env python
# coding: utf-8

# # Package for generation of Ligand-Protein Interaction Fingerprints from MD trajectories (dcd format):
#     IFP include : 
#         PL H-bonds
#         PL salt bridges
#         PL hydrophobic contacs
#         PL interactions that involve aromatic rings
#         PL interactions that involve halogen atom
#         P-water-L bonds
# 
#############################
### v 1.0
#
#    Copyright (c) 2020
#    Released under the EUPL Licence, v1.2 or any higher version
#    
### Author: Daria Kokh
#    Daria.Kokh@h-its.org
#    Heidelberg Institute of Theoretical Studies (HITS, www.h-its.org)
#    Schloss-Wolfsbrunnenweg 35
#    69118 Heidelberg, Germany
################################# 
#
# 
# ### Input data required:
#     - a reference pdb file of the system (for example, generated from the first frame)
#     - ligand mol2 and pdb files
#     
#     
# ### Packages required:
    #     numpy
    #     pandas
    #     matplotlib
    #     seaborn
    #     MDAnalysis        
    #     RDkit
    #     scipy  and  sklearn  kmodes
    #     code is written on Python 3 and tested on the version 3.7
#
# ### Package Overview:
# 
#     IFP_list(property_list, sel_ligands,RE=True) - generation a list of IFP types based on the ligand atom type
#     make_IFT_table(IFP_prop_list,snaps)  - generation of IFP data table 
#     IFP(u_mem,sel_ligands,property_list,WB_analysis = True,RE = True) - generation of an IFP lists
#     table_combine (df_HB,df_WB,df_prop,ligand_name,residues_name = [],start=0,stop=None,step=1)
#     read_IFP(list_IFP)
#     Plot_IFP(df,contact_collection=None)
# 
# 



import glob, os
import sys
import builtins

import numpy as np
import pandas as pd
from pandas import ExcelFile 

from matplotlib import *
from matplotlib import cm
from matplotlib import gridspec
import  pylab as plt
import seaborn as sns

from scipy import stats

from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig

import MDAnalysis as mda
from MDAnalysis.analysis import contacts,align,rms
from MDAnalysis.analysis.base import AnalysisFromFunction
from MDAnalysis.coordinates.memory import MemoryReader
import MDAnalysis.analysis.hbonds as hb
from MDAnalysis.lib.distances import capped_distance, calc_angles

import warnings
warnings.filterwarnings("ignore")

#mda.core.flags['use_periodic_selections'] = True
#mda.core.flags['use_KDTree_routines'] = False


#################################at_negative####
# thresholds:
#####################################
r_cat = 5 # cation-aromatic
r_ari = 5.5  # pi-pi
r_hyd = 4.0 # hydrophobic
r_sar = 4.5 # S-aromatic
r_sal = 4.5 # salt bridge
r_hal = 3.5 # halogen interactions
r_wat = 3.5 # water shell
r_dis = 5.0 # all protein-ligand contacts  
r_lip = 5.0 # specific residues (in particular, lipids)
r_ion = 3.4  # salt bridges with ions

at_aromatic = "((resname PHE TRP TYR HIS HIE HID HE2) and (name CZ* CD* CE* CG* CH* NE* ND*))"
at_positive =  "((resname ARG LYS ) and (name NH* NZ)) or ((resname HI2 ) and (name HD HE))"
at_negative = " ((resname ASP GLU) and (name OE* OD*))"
at_sulfur = "(protein and (name S*))"
at_hydrophob = " (protein  and (name C*  S*) and (not  (name CG and resname ASN ASP))   and (not  (name CD and resname GLU GLN ARG))  and (not  (name CZ and resname TYR ARG))  and (not  (name CE and resname LYS)) and (not  (name CB and resname SER THR))   and (not backbone))"
at_pos_ions = "(resname  MN ZN Mn Zn Ca CA NA Na)"
at_water = "(resname WAT HOH SOL TIP3)"
at_halogens = "( type I CL BR Br Cl)"
at_noH = "( not name H* )"

angle_CHal_O = 150  # currently is not used

resi_aromatic = ["HIS","HIE","HID","HI2","TYR","TRP","PHE"]

#######################################################################
#
#      FUNCTION FOR computing of contact properties for a particular compound and interaction type
#
#######################################################################

class  IFP_prop:
    """
    CLASS of the Interaction Fingerprint properties for particular type of the protein-ligand interection
    
    name - type of interaction
    atoms -ligand atom of that belong to this type
    sel_a - string describing a selection of protein residues/atoms 
    sel_b - string describing a selection of ligand atoms 
    dist - distance to be used for analysis
    contacts - list of contact 
    """
    def __init__(self,name,atoms,sel_a,sel_b,dist):
        self.name = name   
        self.atoms = atoms   
        self.sel_a = sel_a   
        self.sel_b = sel_b    
        self.dist = dist     
        self.contacts = []  


#------------------------

#######################################################################
#
#      FUNCTION taht makes a list for computing specific contacts using MDAnalysis
#
#######################################################################
    
def IFP_list(property_list, sel_ligands, RE=True, Lipids = []):
    """
    Parameters:
    property_list - ligand atom properties as generated by Rdkit
    sel_ligands - ligand residue name
    
    Returns:
    a list of properties that can be then used to extract IFP from a snapshot
    """
    IFP_prop_list = []
    try:  # hydrophobic
        line = ""
        if "Hydrophobe" in property_list.keys():
            for l in tuple(set(property_list["Hydrophobe"])): line = line + l + " "
            sel_a = at_hydrophob
            sel_b = '((resname '+sel_ligands+") and (name "+line+") )"
            IFP_prop_list.append(IFP_prop("HY",line,sel_a,sel_b,r_hyd))
    except:
        print("HY failed")
        pass
    try:  #--- salt bridge with negatively-charged residues 
        line = ""
        if "PosIonizable" in property_list.keys():
            for l in tuple(set(property_list["PosIonizable"])): line = line + l +" "
            sel_a = at_negative #"((resname ASP GLU) and (name OE* OD*)) "
            sel_b = '((resname '+sel_ligands+") and (name "+line+") )"
            IFP_prop_list.append(IFP_prop("IP",line,sel_a,sel_b,r_sal))
    except:
        print("IP failed")
        pass
    
    try: #--- salt bridges with posetively charged residues  
        line = ""
        if "NegIonizable" in property_list.keys():
            for l in tuple(set(property_list["NegIonizable"])): line = line + l +" "
            sel_a = at_positive # "((resname ARG LYS ) and (name NH* NZ)) or ((resname HI2 ) and (name HD HE))"
            sel_b = '((resname '+sel_ligands+") and (name "+line+") )"
            IFP_prop_list.append(IFP_prop("IN",line,sel_a,sel_b,r_sal))
          #--- salt bridges with positively charged ions 
            sel_a = at_pos_ions 
            sel_b = '((resname '+sel_ligands+") and (name "+line+") )"
            IFP_prop_list.append(IFP_prop("IO",line,sel_a,sel_b,r_ion))
    except:
        print("IN failed")
        pass
    
    try:  #--- cation-aril interactions  and aromatic stacking
        """

        line = ""
        if "PosIonizable" in property_list.keys():
            for l in tuple(set(property_list["PosIonizable"])): line = line + l +" "
        if "Aromatic" in property_list.keys():  # pi-pi
            for l in np.asarray(property_list["Aromatic"]): line = line + l +" "
        sel_b = '((resname '+sel_ligands+" ) and (name "+line+"))"
        sel_a = at_aromatic    #"((resname PHE TRP TYR HI2 HIS HIE HID) and (name CZ* CD* CE* CG* CH* NE* ND*)) "   
        if("PosIonizable" in property_list.keys()): 
            IFP_prop_list.append(IFP_prop("AR",line,sel_a,sel_b,r_cat))
        if("Aromatic" in property_list.keys()):
            IFP_prop_list.append(IFP_prop("AR",line,sel_a,sel_b,r_ari))
        """
        line_p = ""
        line_a = ""
        sel_a = at_aromatic    #"((resname PHE TRP TYR HI2 HIS HIE HID) and (name CZ* CD* CE* CG* CH* NE* ND*)) "   
        if "PosIonizable" in property_list.keys():
            for l in tuple(set(property_list["PosIonizable"])): line_p = line_p + l +" "
            sel_b = '((resname '+sel_ligands+" ) and (name "+line_p+"))"
            IFP_prop_list.append(IFP_prop("AR",line_p,sel_a,sel_b,r_cat))
        if "Aromatic" in property_list.keys():  # pi-pi
            for l in np.asarray(property_list["Aromatic"]): line_a = line_a + l +" "
            sel_b = '((resname '+sel_ligands+" ) and (name "+line_a+"))"
            IFP_prop_list.append(IFP_prop("AR",line_a,sel_a,sel_b,r_ari))
        
    except:
        print("AR1 failed")
        pass
    """
    try: #--- S -aromatic 
        sel_b = '((resname '+sel_ligands+") and (type S ))"
        sel_a = at_aromatic #"((resname PHE TRP TYR HI2 HIS HIE HID) and (name CZ* CD* CE* CG* CH* NE* ND*)) "
        IFP_prop_list.append(IFP_prop("AR",line,sel_a,sel_b,r_sar))
    except:
        print("AR2 failed")
        pass
    """
    try: #--- aromatic ligand - cation , amide
        line = ""
        if "Aromatic" in property_list.keys():
            for l in np.asarray(property_list["Aromatic"]): line = line + l +" "
            sel_b = '((resname '+sel_ligands+" ) and (name "+line+") )"
            sel_a = at_positive # +" or "+at_sulfur #"((resname ARG LYS ) and (name NH* NZ*)) or (backbone and name H)  
            IFP_prop_list.append(IFP_prop("AR",line,sel_a,sel_b,r_cat))
    except:
        print("AR3 failed")
        pass
    
    try: #--- halogen bonds with atromatic, negatively charged or backbone carbonyl oxygen
        sel_b = '((resname '+sel_ligands+" ) and "+at_halogens+" )"
        sel_a = at_aromatic +" or "+at_negative+" or "+ " (backbone and name O) "+" or "+at_sulfur #((resname CYS MET) and (type S))"
        IFP_prop_list.append(IFP_prop("HL","HL",sel_a,sel_b,r_hal))
    except:
        print("HL failed")
        pass
    try: # water shell  
#        sel_a = "((sphzone 12.0 resname "+sel_ligands+") and (resname WAT and type O))"
        sel_a = " (resname WAT HOH SOL TIP3 and name O*) "
        sel_b = '(resname '+sel_ligands+" ) and "+at_noH
        IFP_prop_list.append(IFP_prop("WA","HE",sel_a,sel_b,r_wat))
    except:
        print("WA failed")
        pass
    
    if RE:
        try: # any protein-ligand contacts 
            sel_a = "(not resname WAT HOH SOL TIP3) and (not name H*)"
            sel_b = '(resname '+sel_ligands+" ) and "+at_noH
            IFP_prop_list.append(IFP_prop("RE","HE",sel_a,sel_b,r_dis))
        except:
            print("RE failed")
            pass
        
    if len(Lipids) > 0:
        try: # any lipid-ligand contacts 
            line = ""
            for l in Lipids: line = line + l +" "
            sel_a = "((resname '+line+' ) and "+at_noH +") "
            sel_b = '(resname '+sel_ligands+") and "+at_noH
            IFP_prop_list.append(IFP_prop("LL",line,sel_a,sel_b,r_lip))
        except:
            pass
    return (IFP_prop_list)


#######################################################################
#
#      FUNCTION to combine a list of IFP of a particular trajectory in one matrix
#
#######################################################################
def make_IFT_table(IFP_prop_list,snaps,columns_extended = []):
    """
    Most of interections are taken from the list given in 
    https://www.cambridgemedchemconsulting.com/resources/molecular_interactions.html
    
    Parameters:
    IFP_prop_list list of IFP objects 
    
    Returns:
    column names and matrix with IFP    values
    """
 
    # first make a list of all types of contacts observed (columns)
    if len(columns_extended) == 0:
        columns = []
        for IFP_type  in IFP_prop_list:  # loop over frames
            for s in IFP_type.contacts:  # loop over different contacts
                for c in s[1]:           # loop over  particular contacts in a particular frame                
                    if(IFP_type.name == "WA"):
                        IFP_element = "WAT"
                    elif(IFP_type.name == "LL"):
                        IFP_element = "LIP"
                    else:# combine contact type with residues type and name
                        IFP_element = c[0]
                    columns.append(IFP_element)
        columns = np.unique(np.asarray(columns).flatten())
    else:
        columns = columns_extended

    
    times = np.linspace(0,snaps,snaps)
    IFP_matrix=np.zeros((len(times),len(columns)),dtype = np.int8)
    
    for IFP_type  in IFP_prop_list: # loop over different contacts
        for s in IFP_type.contacts: # loop over frames
            for c in s[1]:          # loop over  particular contacts in a particular frame (frame - s[0])
                if(IFP_type.name == "WA"):
                    IFP_element = "WAT"
                if(IFP_type.name == "LL"):
                    IFP_element = "LIP"
                else:# combine contact type with residues type and name
 #                   IFP_element = IFP_type.name+"_"+c[0]+str(c[1]) 
                    IFP_element = c[0]
                try:
                    col = np.argwhere(columns == IFP_element).flatten()[0]
                except:
                    print("ERROR: ",IFP_element," was not found in ",columns,len(columns_extended))
                try:
                    if((IFP_type.name == "WA") ): IFP_matrix[s[0],col] += 1
                    else: IFP_matrix[s[0],col] = 1
                except:
                    print("IFP was not found: ",IFP_element,col,s[0])
    return(columns,IFP_matrix)

    
 #######################################################################
#
#     FUNCTION FOR generation of interaction fingerprints (IFP) in a trajectory
#
#######################################################################
def IFP(u_mem,sel_ligands,property_list, WB_analysis = True, RE = True,Lipids = [],WB_debug = False):
    import datetime
    """
    Parameters:
    u - trajectory - universe object
    ligand name -  ligand residue name
    property_list - python dictionary of ligand atom properties (created by ligand_analysis)
    
    Reterns:
    """

    #---------------------------------------------------------------
    #- find hydrogen bonds between ptotein and ligand  
    donor_list = []
    
    hb.HydrogenBondAnalysis.DEFAULT_DONORS['OtherFF'] = hb.HydrogenBondAnalysis.DEFAULT_DONORS['CHARMM27']
    hb.HydrogenBondAnalysis.DEFAULT_ACCEPTORS['OtherFF'] = hb.HydrogenBondAnalysis.DEFAULT_ACCEPTORS['CHARMM27']
    hb.WaterBridgeAnalysis.DEFAULT_DONORS['OtherFF'] = hb.WaterBridgeAnalysis.DEFAULT_DONORS['CHARMM27']
    hb.WaterBridgeAnalysis.DEFAULT_ACCEPTORS['OtherFF'] = hb.WaterBridgeAnalysis.DEFAULT_ACCEPTORS['CHARMM27']
    if "Donor" in set(property_list):
        donor_line = tuple(set(property_list["Donor"]))
        hb.HydrogenBondAnalysis.DEFAULT_DONORS['OtherFF'] = hb.HydrogenBondAnalysis.DEFAULT_DONORS['OtherFF']+donor_line
        hb.WaterBridgeAnalysis.DEFAULT_DONORS['OtherFF'] = hb.WaterBridgeAnalysis.DEFAULT_DONORS['OtherFF']+donor_line
    if "Acceptor" in set(property_list):      
        acceptor_line = tuple(set(property_list["Acceptor"]))
        hb.HydrogenBondAnalysis.DEFAULT_ACCEPTORS['OtherFF'] = hb.HydrogenBondAnalysis.DEFAULT_ACCEPTORS['OtherFF']+acceptor_line
        hb.WaterBridgeAnalysis.DEFAULT_ACCEPTORS['OtherFF'] = hb.WaterBridgeAnalysis.DEFAULT_ACCEPTORS['OtherFF']+acceptor_line
    
    h = hb.HydrogenBondAnalysis(u_mem, selection1 ='resname '+sel_ligands,selection2=' not resname WAT HOH SOL '+sel_ligands, distance=3.3, angle=100, forcefield='OtherFF')
    print("Start HB analysis",datetime.datetime.now().time())
    h.run()
    h.generate_table()
    df_HB = pd.DataFrame.from_records(h.table)
    
    #---------------------------------------------------------------
    #------ find water bridges between ligand and protein--------
    df_WB = pd.DataFrame()
    if WB_analysis :
        print("Start WB analysis",datetime.datetime.now().time())
#        try:
            # will add O atom of water as an HB donor (it is concidered as acceptor by default assuming as a protein backbone atom)
#            hb.WaterBridgeAnalysis.DEFAULT_DONORS['OtherFF'] = hb.WaterBridgeAnalysis.DEFAULT_DONORS['OtherFF']+tuple(set("O"))        
#            w = hb.WaterBridgeAnalysis(u_mem, selection1  = 'resname '+sel_ligands, selection2  = ' not resname WAT HOH SOL '+sel_ligands,water_selection=" resname WAT HOH SOL ",  distance=3.0, angle=120, forcefield='OtherFF',order=WB_order)
#            w.run()
#            w.generate_table()
#            df_WB = pd.DataFrame.from_records(w.table)
            ################ test
        df_WB = Water_bridges(u_mem, sel_ligands,WB_debug)
            ############### test ele1_index	sele2_index	sele1_resnm	sele1_resid	sele1_atom	sele2_resnm	sele2_resid	sele2_atom
#        except:
#            print("Water bridge analysis failed ")
        
    #---------------------------------------------------------------
    #------ find all other contacts--------   
    print("Start collecting IFPs: ",datetime.datetime.now().time())
    IFP_prop_list = IFP_list(property_list, sel_ligands,RE,Lipids)
    u_list_all = []
    for IFP_type  in IFP_prop_list:
        line = IFP_type.sel_a +" and around "+str(IFP_type.dist)+" "+ IFP_type.sel_b
        try:
            u_list_all.append(u_mem.select_atoms(line, updating=True))
        except:
            print("Selection Error for the type: "+IFP_type.name+" ;  "+line)
    #--------------------------------------------------
    start = 0
    IFPs_unique_list = []
    for i in range(len(u_mem.trajectory)):
        u_mem.trajectory[i] 
        for u_list,IFP_type  in zip(u_list_all,IFP_prop_list):
            found = []
            if (IFP_type.name == "WA"): 
                for u in u_list:   found.append(["WAT",u.name])
            elif (IFP_type.name == "LL"): 
                for u in u_list:   found.append(["LIP",u.name])
            elif (IFP_type.name == "AR"):
                u_ar = []
                u_ar_n = []
                for u in u_list:  
                    u_ar.append(u.resid)
                    u_ar_n.append(u.resname)
                if len(u_ar)> 0:
                    ar_resid, ar_n = np.unique(u_ar,return_counts=True)
                    for u in u_list:  
                        # check if this aromatic residue has more than 4 contacts with an aromatic fragment of a ligand
   #                     print("AROMATIC  ",ar_resid,ar_n,u_ar,u_ar_n)
                        if(u.resid in ar_resid[ar_n > 4]): 
                            # check also residue name to deal the case of residues with the same id
                            if( np.unique(np.asarray(u_ar_n)[np.where(u_ar==u.resid)[0]]).shape[0])== 1: 
                                found.append([IFP_type.name+"_"+u.resname+str(u.resid),u.name])
                        # here we will check if cation (LYS or ARG) really contact an aromatic ring of the ligand
                        elif(u.resname in ["LYS","ARG"]):
                            cation = u.resname 
                            if("Aromatic" in property_list.keys()):
                                line_ar = ""
                                for l in np.asarray(property_list["Aromatic"]): line_ar = line_ar + l +" "
                                line1 = "(resname "+sel_ligands+" and ( not type H O) and name "+line_ar+") and around "+str(r_cat)+" (resid "+str(u.resid[u.resname == cation][0]) + " and type N )" 
                                u1_list = (u_mem.select_atoms(line1,updating=True))
                                if(len(u1_list) > 4): 
                                    found.append([IFP_type.name+"_"+u.resname+str(u.resid),u.name])
    #                            print("!!!!!!!!!!!!!  Cat-Ar  interactions found",len(u1_list),u1_list)
                        # now we check if aromatc residue is perpendicular to the aromatic fragment of the ligand (face-to-edge)
                        ### TOBE checked if this works!!!!!================================================
                        elif(u.resname in resi_aromatic) and (u.resid in ar_resid[ar_n <= 4]):
                            if("Aromatic" in property_list.keys()):
                                line_ar = ""
                                for l in np.asarray(property_list["Aromatic"]): line_ar = line_ar + l +" "
                                line1 = "(resname "+sel_ligands+" and ( not type H O) and name "+line_ar+") and around "+str(r_ari)+" (resid "+str(u.resid) + " and (name NE* ND* CE* CD* CZ* CH* CG*))" 
                                u1_list = (u_mem.select_atoms(line1,updating=True))
                                if(len(u1_list) > 4): 
                                    found.append([IFP_type.name+"_"+u.resname+str(u.resid),u.name])
                        ### TOBE checked if this works!!!!!================================================
            elif (IFP_type.name == "HL"):
                for u in u_list:   
                    # here we will check if the angle  C-HAL.... O is about 180grad fro this we look what atoms within r_hal+1 from O - should be only Hal
                    if (u.type == "O") or (u.type == "S"):
                        line1 ="(resname "+sel_ligands+" ) and around "+str(r_hal+1.0)+" (protein and resid "+str(u.resid)+" and name O* S* )"
                        u1_list = (u_mem.select_atoms(line1,updating=True))
#                        print(u.resid,u.name,":::",len(u1_list),u1_list)
                        if len(u1_list) < 2:
                            found.append([IFP_type.name+"_"+u.resname+str(u.resid),u.name]) 
 #                       else: print("HL-O contact found but will not be counted because of the too small angle: ",len(u1_list),u.resname+str(u.resid),u.name)
                        # TO BE DONE instead of previous criterion
                        """
                        else:
                            u1_list = (u_mem.select_atoms(" (resid "+str(u.resid)+" and type O )",updating=True)) # resi
                            u2_list = (u_mem.select_atoms("(resname "+sel_ligands+" and type Cl CL Br BR I) and around "+str(r_hal+1.0)+" (resid "+str(u.resid)+" and type O)",updating=True)) # ligand hal
                            u3_list = (u_mem.select_atoms("(resname "+sel_ligands+" and type C) and around "+str(r_hal+1.0)+" (resid "+str(u.resid)+" and type O)",updating=True)) # ligand carbon
                            B_center = B.centroid(u1_list)
                            BA = A.centroid(u2_list) - B_center
                            BC = C.centroid(u3_list) - B_center
                            alpha = np.arccos(np.dot(BA, BC)/(norm(BA)*norm(BC)))
                            if alpha > angle_CHal_O:
                                found.append([IFP_type.name+"_"+u.resname+str(u.resid),u.name])
                        """

                # now we check if halogen atom is  perpendicular to the aromatic residue
                 ### TOBE checked if this works!!!!!====HAL=====================================
                u_ar = []
                [u_ar.append(u.resid) for u in u_list if u.resname in resi_aromatic]
                if len(u_ar)> 0:
                    ar_resid, ar_n = np.unique(u_ar,return_counts=True)
                    if(u.resid in ar_resid[ar_n > 4]): 
                            found.append([IFP_type.name+"_"+u.resname+str(u.resid),u.name])  
#                    else: print("HL-aromatic contact found but will not be counted because it has less than 4 contacts with aromatic atoms:",ar_resid,ar_n)
                 ### TOBE checked if this works!!!!!====HAL========================================
            else:  
                #print("HY/IP: ",IFP_type.name,len(u_list))
                for u in u_list:
                        found.append([IFP_type.name+"_"+u.resname+str(u.resid),u.name]) 
#                        if IFP_type.name == "HL": print("HAL:",u.resname,u.resid,u.name)
                
            if(found): 
                IFP_type.contacts.append((i,found))
                if start == 0:  
                    IFPs_unique_list = np.unique(np.asarray(found)[:,0])
                    start += 1
                else:  
                    IFPs_unique_list = np.unique(np.append(IFPs_unique_list,np.asarray(found)[:,0]))
#                print(IFPs_unique_list)

    print("Start building IFP table: ",datetime.datetime.now().time())
    if(len(IFP_prop_list) > 0):
        columns,IFP_matrix = make_IFT_table(IFP_prop_list,len(u_mem.trajectory),columns_extended = IFPs_unique_list)
        df_prop = pd.DataFrame( data = IFP_matrix,index=None, columns=columns) 
    else:
        print("Something is wrong - IFP property list is empty")
    print("IFP database is ready ",datetime.datetime.now().time())
    return(df_prop,df_HB,df_WB) 

####################################################################################################
#
# combine three tables 
#    - with protein-ligand hydrogen bonds
#    - with protein-ligand water bridges
#    - with other protein-ligand interaction properties(IFP)
# in one table
#
####################################################################################################
def table_combine(df_HB,df_WB,df_prop,ligand_name,residues_name = [],start=0,stop=None,step=1):
    """
    Parameters:
    df_HB - H-bond table
    df_prop - IFP table 
    ligand_name      - ligand nale
    residues_name a list of properties that (column names) that can be used to generate tables with the same column list
    
    Return:
    updated table
    """
    if stop :
        if len(range(start,stop,step)) != np.asarray(df_prop.shape[0]):
            stop = (df_prop.shape[0] -start)*step
    else:
        stop = df_prop.shape[0]
        
#---------------- extract hydrogen bonds between ligand and protein and add to IFP table----------        
    columns_resname = []
    df_prop["time"] =range(start,stop,step)
    #------- list of the residues making donor-HB  with the  ligand, but not water
    df_noWatD = df_HB[~df_HB.donor_resnm.isin([ligand_name,"WAT"])]  # protein donor
    df_noWatD = df_noWatD[df_noWatD.acceptor_resnm == ligand_name]   # protein donor and ligand acceprot
    
    #------- list of the residues making acceptor-HB  with the  ligand , but not water
    df_noWatA = df_HB[~df_HB.acceptor_resnm.isin([ligand_name,"WAT"])]  # protein acceptor
    df_noWatA = df_noWatA[df_noWatA.donor_resnm == ligand_name]  # protein acceptor and ligand donor

    t_list = []    
    for t in df_HB.time.unique().tolist():
        raw = int(t)
        if not df_noWatD.empty:
            df_noWatD_t = df_noWatD[(df_noWatD.time == t)]
            if not df_noWatD_t.empty:
                for d in df_noWatD_t.donor_resid.tolist():
                    r = "HD_"+df_noWatD_t[df_noWatD_t.donor_resid == d].donor_resnm.tolist()[0]+str(d)
                    if r not in columns_resname:  columns_resname.append(r)
                    t_list.append((raw,r))
        if not df_noWatA.empty:
            df_noWatA_t = df_noWatA[(df_noWatA.time == t)]
            if not df_noWatA_t.empty:
                for d in df_noWatA_t.acceptor_resid.tolist():
                    r = "HA_"+df_noWatA_t[df_noWatA_t.acceptor_resid == d].acceptor_resnm.tolist()[0]+str(d)
                    if r not in columns_resname:  columns_resname.append(r)
                    t_list.append((raw,r))
    properties =  np.zeros((len(df_prop.index.values.tolist()),len(columns_resname)),dtype=np.int8)
    for j,c in enumerate(np.sort(np.asarray(columns_resname))):
        for i,cc in enumerate(t_list):
            if(c == cc[1]):   
                properties[cc[0],j] = 1
        df_prop[c] =   properties[:,j]    
#    print("--------------HB found------\n")
    
            
#---------------- extract water bridges between ligand and protein and add to IFP table----------
    if not df_WB.empty:
        # we have to check naming since it was changed in later version
        new_df_WB_columns = []
        info = True
        #------ new version 16-12-2019
        t_list = []
        column_resi = []
        # get a list of INH-WAT (sele1 - sele2)
        df_WB_INH = df_WB[(df_WB.sele1_resnm.isin([ligand_name]) & df_WB.sele2_resnm.isin(["WAT","HOH","SOL","TIP3"]))]
#        print(df_WB_INH)
        # get a list of WAT-Prot (sele1 - sele2)
        df_WB_Prot = df_WB[(~(df_WB.sele2_resnm.isin([ligand_name,"WAT","HOH","SOL","TIP3"])) & (df_WB.sele1_resnm.isin(["WAT","HOH","SOL","TIP3"])))]
#        print(df_WB_Prot)
 #       df_WB_WAT = df_WB[((df_WB.sele2_resnm == "WAT") & (df_WB.sele1_resnm == "WAT"))]
        for t in df_WB.time.unique().tolist():
            raw = int(t)
            df_WB_t = df_WB[df_WB.time == t]
          #  df_WB_WAT_t = df_WB_WAT[df_WB_WAT.time == t]
            df_WB_Prot_t = df_WB_Prot[df_WB_Prot.time == t]
            df_WB_INH_t = df_WB_INH[df_WB_INH.time == t]
            if ((not df_WB_Prot_t.empty) and (not df_WB_INH_t.empty)):
                common_WAT = np.intersect1d(df_WB_INH_t.sele2_resid.values,df_WB_Prot_t.sele1_resid.values)
                for r in np.unique(df_WB_Prot_t.sele2_resid.values):
                    r1 = "WB_"+df_WB_Prot_t[df_WB_Prot_t.sele2_resid == r].sele2_resnm.values[0]+str(r)
                    t_list.append((raw,r1))
                    if r1 not in column_resi: 
                        column_resi.append(r1)
          #     print(int(t),"---> common WAT --->",common_WAT,"Prot: ",np.unique(df_WB_Prot_t.sele2_resid.values))
        #------------------------------
        properties =  np.zeros((len(df_prop.index.values.tolist()),len(column_resi)),dtype=np.int8)
        for j,c in enumerate(column_resi):
            for i,cc in enumerate(t_list):
                if(c == cc[1]):   
                    properties[cc[0],j] = 1
            df_prop[c] =   properties[:,j]    
    
#-----------------------------------------------------------------    
    # add more columns for IFP provided as input but not found in the current--
    for rr in residues_name:
        if (rr in df_prop.columns.tolist()):   pass
        else:  df_prop[rr] = 0
            

#----------------------------------------------------------------        
    # cleaning the table 
    # fist order residues by number
    df_prop_order_new = []
    for df_prop_order in df_prop.columns.tolist():
        if df_prop_order.find("_") > 0:
            df_prop_order_new.append(int(df_prop_order[df_prop_order.find("_")+4:]))
        else: df_prop_order_new.append(0)            
    properties = np.asarray(df_prop.columns.tolist())[np.argsort(df_prop_order_new)]
    # then use order of the properties HY - HD - HA - IP - IN  
    for i_df in range(1,len(properties)):
        if (properties[i_df].find("_")> 0) and (properties[i_df-1].find("_")> 0):
            if(properties[i_df][properties[i_df].find("_"):] == properties[i_df-1][properties[i_df-1].find("_"):]):
                properties[i_df-1:i_df+1] = np.sort(properties[i_df-1:i_df+1])
#        print("---",properties)
    df_prop = df_prop[properties]
    #--- change column position puttinhg time at the beginning and Water at the end---
    df_prop = df_prop[np.concatenate((["time"],df_prop.columns[df_prop.columns != "time"].tolist()))]
    if "WAT" in df_prop.columns.tolist():
        df_prop = df_prop[np.concatenate((df_prop.columns[df_prop.columns != "WAT"].tolist(),["WAT"]))]
    
    return (df_prop)

#########################################################################
# Program for reading IFP databases
# Additionally column with ligand name is added
# and COM column is splitted to COM_x, COM_y, COM_z
########################################################################
def read_IFP(list_IFP):
    """
    Parameters:
    dictionary of files with ITP databases {name1:file_path1[,name2:filepath2],...}
    
    Returns:
    combined IFP database
    """
    unpickled_dfi = []
    ligandsi = []
    for lig in list_IFP:
        unpickled_dfi.append(pd.read_pickle(list_IFP[lig]))
        ligandsi.append(lig)

    # add ligand names and make a joint list of columns
    intersect = []
    for (df,lig) in zip(unpickled_dfi,ligandsi):
        df["ligand"] = np.repeat(lig,df.shape[0]) 
        diff = np.setdiff1d(np.asarray(df.columns.tolist()),intersect)
        if(len(intersect) == 0):  intersect = diff
        else: intersect = np.append(intersect,diff)
    
    # add empty columns for those that are present in the joint list but absent in the database
    unpickled_df = pd.DataFrame(columns=intersect) 
    for (df,lig) in zip(unpickled_dfi,ligandsi):
        for ifp in intersect:
            if ifp not in df.columns.tolist():
                df[ifp] = np.repeat(np.int8(0),df.shape[0])    
        unpickled_df = pd.concat([unpickled_df, df], axis=0, sort=False)

    # converge COM string to  x y z components
    if "COM"  in unpickled_df.columns.tolist():
        COM_x = []
        COM_y = []
        COM_z = []
        for l in unpickled_df.COM:
            COM_x.append(l[0])
            COM_y.append(l[1])
            COM_z.append(l[2])
        unpickled_df["COM_x"] = COM_x
        unpickled_df["COM_y"] = COM_y
        unpickled_df["COM_z"] = COM_z
    
    return(unpickled_df)


########################################
#
#     PLOT IFP  and ligand water shell for a trajectory
#
########################################
def Plot_IFP(df,contact_collection=None,out_name="",ifp_list = ["HY","AR","HD","HA","HL","IP","IN","WB","LIP"]):
    """
    Parameters:
    df- IFP database
    contact_collection - set of IFP to be shown
    
    Returns:
    """
 #   color = ['r','b','forestgreen','lime','m','c','teal','orange','yellow','goldenrod','olive','tomato','salmon','seagreen']
    top = cm.get_cmap('Oranges_r', 128)
    bottom = cm.get_cmap('Blues', 128)
    color = np.vstack((top(np.linspace(0, 1, 128)),
                       bottom(np.linspace(0, 1, 128))))
    columns_IFP = []  # standard IFP
    columns_CONT = []  # just contacts
    for c in df.columns.tolist():
        if c[0:2] in ifp_list:
            columns_IFP.append(c)
        elif c[0:2]  == "RE":
            columns_CONT.append(c)
                
    df1 = df[columns_IFP].values
    if df1.shape[0] < 2: return
    fig = plt.figure(figsize=(16, 6))
    gs = gridspec.GridSpec(1, 3, width_ratios=[4, 2, 1]) 
    ax = plt.subplot(gs[0])
   
    ax =  plt.subplot(gs[0])
    ax.set_title('IFP')
    xticklabels = columns_IFP
    if len(columns_IFP) < 25:  sns.heatmap(np.float32(df1), cmap="YlGnBu", xticklabels=xticklabels)
    else:  sns.heatmap(np.float32(df1), cmap="YlGnBu")
    ax.set_title('IFP')
    
    if( df[columns_CONT].shape[1] > 0):
        ax = plt.subplot(gs[1])
        ax.set_title('RE')
        df1 = df[columns_CONT].values
        sns.heatmap(np.float32(df1), cmap="YlGnBu")
        ax.set_title('Contacts')
    
    if("WAT"  in df.columns.tolist()):
        ax = plt.subplot(gs[2])
        ax.set_ylim(0,max(df["WAT"].tolist()))
        if ("Repl"  in df.columns.tolist()):
            for i,r in  enumerate(np.unique(df.Repl.tolist())):
                plt.plot(df[df.Repl == r]["WAT"], marker='o', linewidth = 0,color = color[i],label=r)
            if np.unique(df.Repl.tolist()).shape[0] < 10:
                plt.legend()
        else:
            plt.plot(df["WAT"],'go')
        ax.set_title('Water shell')
        ax.set_xlabel('frame')
        ax.set_ylabel('# of water molecules')

    if out_name == "" and builtins.SHOULD_SHOW_PLOTS:
        plt.show()
    else:
        plt.savefig(out_name,dpi=300)

    return


###################################
#
##################################
def rank_IFP_resi(df,ifp_type=['AR','HY','HA','HD','HL','IP','IN',"IO","WB"]):
    """    
    script extracts and ranks by the residue number IFP list  from the IFP table 
    
    Parameters:
    df - pkl table
    ifp_type - list of IFP types to be considered
    
    Results:
    columns_IFP - list of IFP
    columns_R
    """
    columns_IFP = []  # standard IFP
    number = []
    for c in df.columns.tolist():
        if c[0:2] in ifp_type: 
            columns_IFP.append(c)
            if c[3:].isdigit():    number.append(int(c[3:]))
            elif c[4:].isdigit():    number.append(int(c[4:]))
            elif c[5:].isdigit():    number.append(int(c[5:]))
            else: number.append(int(c[6:]))
    columns_IFP = np.asarray(columns_IFP)[np.argsort(np.asarray(number))]

    columns_RE = []  # standard IFP
    number = []
    for c in df.columns.tolist():
        if c[0:2] == "RE": 
            columns_RE.append(c)
            if c[3:].isdigit():    number.append(int(c[3:]))
            elif c[4:].isdigit():    number.append(int(c[4:]))
            elif c[5:].isdigit():    number.append(int(c[5:]))
            else: number.append(int(c[6:]))

    columns_RE = np.asarray(columns_RE)[np.argsort(np.asarray(number))]
    return(columns_IFP,columns_RE)


###################################
#
##################################

def Plot_IF_trajectory(df_tot,ifp_type = np.asarray(['AR','HA','HD','HY','IP','IN',"IO","WB"]),head_tail=-1,save_file = ""):
    
    columns_IFP,columns_RE= rank_IFP_resi(df_tot, ifp_type)
    columns_sel = columns_IFP
    if head_tail < 0:
        head_tail = int(df_tot.shape[0]/3)       
    X = []
    df = df_tot[columns_sel] 
    columns_sel = columns_sel[df.mean().values > 0.01]
    columns_HB=[]
    [columns_HB.append(c) for c in columns_sel if (c[:2] == "HD" or c[:2] == "HA")]
    columns_HB = np.asarray(columns_HB)    

    fig = plt.figure(figsize=(14,3),dpi=150)

    df = df_tot[np.append(columns_sel,"time")] 
    n_hb = len(columns_HB[(df[columns_HB].mean()> 0.75).values])
    plt.bar(range(0,len(columns_sel)),df[df.time < head_tail][columns_sel].mean(),alpha=0.6,label="first "+str(head_tail)+" frames")
    plt.bar(range(0,len(columns_sel)),df[df.time > df.shape[0]-head_tail][columns_sel].mean(),alpha=0.6,label="last "+str(head_tail)+" frames")
    plt.bar(np.asarray(range(0,len(columns_sel))),df[columns_sel].mean(),color="",label="all frames",edgecolor ='k',hatch="/")
    plt.xticks(range(0,len(columns_sel)), columns_sel, rotation='vertical',fontsize=10)
    plt.legend(fontsize=10, loc = 'upper left')
    
    if save_file != "":
        plt.savefig(file_save,format='png', dpi=300, bbox_inches='tight',transparent=True)

    if builtins.SHOULD_SHOW_PLOTS:
        plt.show()

    return

#############################################
#   
#
############################################

def Water_bridges(u_mem, sel_ligands, WB_debug = False):
    """
    A very simple procedure for detection of possible protein-ligand water bridges 
    Parameters:
    u_mem - trajectory
    residues_name a list of properties that (column names) that can be used to generate tables with the same column list
    sel_ligands - ligand residue name 
    Returns:
    df_WB - pkl table with all components of water bridges
    """
    col_exchange = {"sele1_index": "sele2_index", "sele2_index": "sele1_index" ,"sele1_resnm": "sele2_resnm", "sele1_resid": "sele2_resid", "sele1_atom": "sele2_atom","sele2_resnm": "sele1_resnm", "sele2_resid": "sele1_resid","sele2_atom": "sele1_atom"}
    col_transfer = {"donor_index": "sele1_index", "acceptor_index": "sele2_index","donor_resnm": "sele1_resnm", "donor_resid": "sele1_resid","donor_atom": "sele1_atom", "acceptor_resnm": "sele2_resnm","acceptor_resid": "sele2_resid","acceptor_atom": "sele2_atom"}
        
    angle_th =100
    dist_th =3.3
    
    #-------------------------------------------------------------
    def clean_dataset(wb_check,sel_ligands):
        """
        function that checks if water molecule has contacts to both protein and ligand and removes water that are not
        Parameters:
        dataset of potential water bridges
        Returns:
        clean dataset
        """  
#        print(wb_check)
        water_list_d = wb_check[wb_check["donor_resnm"].isin(["WAT", "HOH", "SOL","TIP3"])].donor_resid.values
        water_list_a = wb_check[wb_check["acceptor_resnm"].isin(["WAT", "HOH", "SOL","TIP3"])].acceptor_resid.values
        l = np.concatenate((water_list_d, water_list_a))
        res_w, c_w = np.unique(l,return_counts=True)
        fexcl = res_w[c_w<2]
        wb_check_cleaned = wb_check
        # check if water appeares just once - remove it then:
        if fexcl.shape[0] > 0:
            if not wb_check_cleaned[wb_check_cleaned.donor_resid.isin(fexcl)].empty:
                wb_check_cleaned = wb_check_cleaned[~wb_check_cleaned.donor_resid.isin(fexcl)]
            if not wb_check_cleaned[wb_check_cleaned.acceptor_resid.isin(fexcl)].empty:
                wb_check_cleaned = wb_check_cleaned[~wb_check_cleaned.acceptor_resid.isin(fexcl)]
        # check if water connects to different residues:
        for wat in res_w:
            nowat_list_d = wb_check_cleaned[wb_check_cleaned["donor_resid"]== wat]
            nowat_list_a = wb_check_cleaned[wb_check_cleaned["acceptor_resid"] == wat]
            res, c = np.unique(np.concatenate((nowat_list_d.acceptor_resid.values, nowat_list_a.donor_resid.values)),return_counts=True)
            both_prot = True
            #check if ligand has contacts with water
            for r in res: 
                if nowat_list_d[nowat_list_d.acceptor_resid == r].acceptor_resnm.values.shape[0]> 0:
                    if sel_ligands in nowat_list_d[nowat_list_d.acceptor_resid == r].acceptor_resnm.values:  both_prot = False
                if nowat_list_a[nowat_list_a.donor_resid == r].donor_resnm.values.shape[0]> 0:
                    if sel_ligands in nowat_list_a[nowat_list_a.donor_resid == r].donor_resnm.values: both_prot = False
            if res.shape[0] < 2  or both_prot:
                if not wb_check_cleaned[wb_check_cleaned["donor_resid"]== wat].empty:
                    wb_check_cleaned = wb_check_cleaned[~(wb_check_cleaned["donor_resid"]== wat)]
                if not wb_check_cleaned[wb_check_cleaned["acceptor_resid"] == wat].empty:
                    wb_check_cleaned = wb_check_cleaned[~(wb_check_cleaned["acceptor_resid"] == wat)]
#        print(wb_check_cleaned)
        return(wb_check_cleaned)
        #-------------------------------------------------
    
    
    df_WB = pd.DataFrame()
    
    # 1. we will find a list of ligand- water h-bonds in the all trajectory
    hb.HydrogenBondAnalysis.DEFAULT_DONORS['OtherFF'] = hb.WaterBridgeAnalysis.DEFAULT_DONORS['OtherFF']+tuple(set("O"))        
    h = hb.HydrogenBondAnalysis(u_mem, selection1 ='resname '+sel_ligands,selection2='  resname WAT HOH SOL TIP3', distance=dist_th, angle=angle_th, forcefield='OtherFF') #,update_selection1= False)
    h.run()
    h.generate_table()
    wb1_tot = pd.DataFrame.from_records(h.table)
    if WB_debug:  print("1, Lig-Wat"," -----------------\n",wb1_tot)
# 2. make a list of water molecules 
    lista = []
    if not wb1_tot[wb1_tot.donor_resnm == sel_ligands].empty:
        lista = wb1_tot[wb1_tot.donor_resnm == sel_ligands].acceptor_resid.values
    listd = []
    if not wb1_tot[wb1_tot.acceptor_resnm == sel_ligands].empty:
        listd = wb1_tot[wb1_tot.acceptor_resnm == sel_ligands].donor_resid.values
    
    if (len(lista)+ len(listd) > 0):    
        list_wb = "resid "
        for l in np.unique(lista): list_wb = list_wb + str(l)+" " 
        for l in np.unique(listd): list_wb = list_wb + str(l)+" " 
#3. make a table of water- protein contacts (only selected water molecules are considered)
        h = hb.HydrogenBondAnalysis(u_mem, selection1 =list_wb,selection2=' not resname WAT HOH SOL TIP3'+sel_ligands, distance=dist_th , angle=angle_th, forcefield='OtherFF')
        h.run()
        h.generate_table()
        wb2_tot = pd.DataFrame.from_records(h.table)
        if WB_debug:  print("2, Prot-Wat"," -----------------\n",wb2_tot)
            
        if wb2_tot.shape[0] > 0: 
            # loop over frames
            for time in range(len(u_mem.trajectory)):
                u_mem.trajectory[time] 
                wb2 = wb2_tot[wb2_tot.time.astype(int) == time]
                wb1 = wb1_tot[wb1_tot.time.astype(int) == time]
                if wb2.shape[0] > 0:
                    # exclude the cases where the same hydrogen is a donor for both protein and ligand
                    wb2 = wb2[~(wb2.donor_index.isin(wb1[wb1["donor_resnm"].isin(["WAT", "HOH", "SOL","TIP3"])].donor_index))]
                    # exclude the cases where the same oxigen is an acceptor for both protein and ligand
                    wb2 = wb2[~(wb2.acceptor_index.isin(wb1[wb1["acceptor_resnm"].isin(["WAT", "HOH", "SOL","TIP3"])].acceptor_index))]
                    wb12 =wb1.append(wb2)
# 5. check additionally angles                    
# 5(a) make a list water molecules  that have H-bonds with a ligand
#                    list_ld = wb12[wb12.acceptor_resnm ==  sel_ligands].donor_resid.values
#                    list_la = wb12[wb12.donor_resnm == sel_ligands].acceptor_resid.values
#                    list_w_l = np.concatenate((list_la, list_ld))
#                    if len(list_w_l) == 0: continue 
#                    wb12 = wb12[(wb12.donor_resid.isin(list_w_l))].append(wb12[(wb12.acceptor_resid.isin(list_w_l))])
                    wb12 = clean_dataset(wb12,sel_ligands)
                    if WB_debug: print(time," -----------------\n",wb12)
                    if wb12.empty: continue
                    wat_donor = wb12[wb12["donor_resnm"].isin(["WAT", "HOH", "SOL","TIP3"])]
                    wat_acceptor = wb12[wb12["acceptor_resnm"].isin(["WAT", "HOH", "SOL","TIP3"])]                   
                    lig_donor = wb12[wb12["donor_resnm"].isin([sel_ligands])]
                    prot_donor = wb12[~(wb12["donor_resnm"].isin(["WAT", "HOH", "SOL","TIP3", sel_ligands]))]
                    
# 5(a) ---- check the angle for prot/lig(H)...water(O)...water(H) and remove lines in the table  where it is less then threshold
                    if (not prot_donor.empty)  or (not lig_donor.empty):
                        check_list =[]
                        # loop over water donor atoms
                        for wat_hid,wat_hat in zip(wat_donor.donor_resid.values,wat_donor.donor_atom.values):
                        #    if wat_hid not in list_w_l: continue
                            wat_oid = wat_acceptor[wat_acceptor["acceptor_resid"] == wat_hid].acceptor_resid.values 
                            # loop over water acceptor atoms
                            for wid in wat_oid:
                                # loop over non-water donor bound to water acceptors
                                wid_nonwat = wat_acceptor[wat_acceptor["acceptor_resid"] == wid]
                                for (nowat_hid,nowat_hat,nowat_hin) in zip(wid_nonwat.donor_resid,wid_nonwat.donor_atom,wid_nonwat.donor_index):
                                        if (wat_hid,wat_hat,nowat_hid,nowat_hat,nowat_hin) not in check_list:
                                            check_list.append((wat_hid,wat_hat,nowat_hid,nowat_hat,nowat_hin))
 #                       print("Check list:  ", check_list)
                        for (wat_hid,wat_hat,nowat_hid,nowat_hat,nowat_hin) in check_list:
                            try:
                                    angles = np.rad2deg(
                                        calc_angles(
                                            u_mem.select_atoms("resid "+str(wat_hid)+" and name "+wat_hat,updating=True).positions,
                                            u_mem.select_atoms("resid "+str(wat_hid)+" and name O*",updating=True).positions,
                                            u_mem.select_atoms("resid "+str(nowat_hid)+" and name "+nowat_hat,updating=True).positions,
                                        )
                                        )
#                                    print("Check angle: ",np.round(angles[0],1)," resid "+str(wat_hid)+" "+wat_hat," resid "+str(wid)+" O","resid "+str(nowat_hid)+" "+nowat_hat)
                            except:
                                    angles = 0     
#                                    print("Warning: problem with WB angles (maybe some residue numbers are duplicated): "+" resid "+str(wat_hid)+" "+wat_hat," resid "+str(wid)+" O","resid "+str(nowat_hid)+" "+nowat_hat)
                            if angles < angle_th:   
                                        wr =(wat_hid,nowat_hin)
                                        if WB_debug:  print("REMOVE incorrect H-bonds from the WB list:",wr,nowat_hid,nowat_hat)
                                        wb12 = wb12[~((wb12.acceptor_resid == wr[0]) &  (wb12.donor_index == wr[1]))]
#                                        print("INTERMEDIATE -----------------\n",wb12[( (wb12.acceptor_resid == wr[0]) & (wb12.donor_index == wr[1]))])
        
# 5(b) make a list water molecules  that have H-bonds with a ligand, but first revise table
#                    print(time," -----------------\n",wb12)
#                    tt = clean_dataset(wb12)
#                    list_ld = wb12[wb12.acceptor_resnm ==  sel_ligands].donor_resid.values
#                    list_la = wb12[wb12.donor_resnm == sel_ligands].acceptor_resid.values
#                   list_w_l = np.concatenate((list_la, list_ld))
#                    wb12 = wb12[(wb12.donor_resid.isin(list_w_l))].append(wb12[(wb12.acceptor_resid.isin(list_w_l))])
#                    if len(list_w_l) == 0: continue        
                    wb12 = clean_dataset(wb12,sel_ligands)
                    if wb12.empty: continue
#                    print(time," -----------------\n",wb12)
                    
                    wat_donor = wb12[wb12["donor_resnm"].isin(["WAT", "HOH", "SOL","TIP3"])]                    
                    lig_acceptor = wb12[wb12["acceptor_resnm"].isin([sel_ligands])]
                    prot_acceptor = wb12[~(wb12["acceptor_resnm"].isin(["WAT", "HOH", "SOL","TIP3", sel_ligands]))]
                    
# 5(b) ---- check the angle for prot/lig(O)...water(H)...water(O) 
                    if (not prot_acceptor.empty)  or (not lig_acceptor.empty):
                        check_list =[]
                        # loop over water acceptor atoms
                        for wat_hid,wat_hat in zip(wat_donor.donor_resid.values,wat_donor.donor_atom.values):
                            oid = wat_donor[wat_donor["donor_resid"] == wat_hid].acceptor_resid.values 
                            oin = wat_donor[wat_donor["donor_resid"] == wat_hid].acceptor_index.values 
                            oat = wat_donor[wat_donor["donor_resid"] == wat_hid].acceptor_atom.values
                            for (oid_nonwat, oat_nonwat, oin_nonwat) in zip(oid,oat,oin):
                                if( wat_hid,wat_hat,oid_nonwat, oat_nonwat, oin_nonwat) not in check_list :
                                    check_list.append((wat_hid,wat_hat,oid_nonwat, oat_nonwat, oin_nonwat))
#                        print("Check list: ",check_list)
                        for (wat_hid,wat_hat,oid_nonwat, oat_nonwat, oin_nonwat) in check_list:
                            try:
                                    angles = np.rad2deg(
                                        calc_angles(
                                            u_mem.select_atoms("resid "+str(oid_nonwat)+" and name "+oat_nonwat,updating=True).positions,
                                            u_mem.select_atoms("resid "+str(wat_hid)+" and name "+wat_hat,updating=True).positions,
                                            u_mem.select_atoms("resid "+str(wat_hid)+" and name O*",updating=True).positions,
                                        )
                                        )
#                                    print("Check angle: ",np.round(angles[0],1)," resid "+str(wat_hid)+" "+wat_hat," resid "+str(wid)+" O","resid "+str(oid_nonwat)+" "+oat_nonwat)
                            except:
                                    angles = 0     
#                                    print("Warning: problem with WB angles (maybe some residue numbers are duplicated): "+" resid "+str(wat_hid)+" "+wat_hat," resid "+str(wid)+" O","resid "+str(oid_nonwat)+" "+oat_nonwat)
                            if angles < angle_th:   
                                        wr =(wat_hid,wat_hat,oin_nonwat)
                                        if WB_debug:  print("REMOVE incorrect H-bonds from the WB list:",wr,oid_nonwat, oat_nonwat)
                                        wb12 = wb12[~((wb12.donor_resid == wr[0]) & (wb12.donor_atom == wr[1]) & (wb12.acceptor_index == wr[2]))]
                            
                                                       
#6. change structure - rename donor/acceptor to sel_1 / sel_2 and move ligand to sel_1                            
                    wb12 =wb12.rename(columns=col_transfer)
                    # bring table to the following structure: sel1 - lig or water ; sel2- protein or water
                    # ---- sel1: lig - sel2: water
                    WB_t = wb12[wb12.sele1_resnm.isin([sel_ligands])] 
                    # ---- sel1: wat - sel2: any
                    wb12_t = wb12[wb12.sele1_resnm.isin(["WAT", "HOH", "SOL","TIP3"])] 
                    # ---- sel1: wat - sel2: prot
                    WB_t = WB_t.append(wb12_t[~(wb12_t.sele2_resnm.isin([sel_ligands]))])  
                    # ---- sel1:water - sel2: lgand  - replace sel1 and sel2
                    WB_t = WB_t.append(wb12[wb12.sele2_resnm.isin([sel_ligands])].rename(columns=col_exchange))
                    # ---- sel1:prot - sel2: wat - replace sel1 and sel2
                    WB_t = WB_t.append(wb12[~wb12.sele1_resnm.isin(["WAT", "HOH", "SOL","TIP3",sel_ligands])].rename(columns=col_exchange)) 
                    df_WB = df_WB.append(WB_t)
 #   print(df_WB)
    return df_WB

