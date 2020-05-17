
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
#    Released under the GNU Public Licence, v2 or any higher version
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


import warnings
warnings.filterwarnings("ignore")

mda.core.flags['use_periodic_selections'] = True
mda.core.flags['use_KDTree_routines'] = False



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
            sel_a = " (protein and (type C  S) and (not  (name CG and resname ASN ASP))   and (not  (name CD and resname GLU GLN ARG))  and (not  (name CZ and resname TYR ARG))  and (not  (name CE and resname LYS)) and (not  (name CB and resname SER THR))   and (not backbone))"                 
            sel_b = '((resname '+sel_ligands+") and (name "+line+") )"
            IFP_prop_list.append(IFP_prop("HY",line,sel_a,sel_b,4.0))
    except:
        print("HY failed")
        pass
    try:  #--- salt bridge with posetively-charged residues 
        line = ""
        if "PosIonizable" in property_list.keys():
            for l in tuple(set(property_list["PosIonizable"])): line = line + l +" "
            sel_a = "((resname ASP GLU) and (name OE* OD*)) "
            sel_b = '((resname '+sel_ligands+") and (name "+line+") )"
            IFP_prop_list.append(IFP_prop("IP",line,sel_a,sel_b,4.5))
    except:
        print("IP failed")
        pass
    
    try: #--- salt bridges with negatively charged residues  
        line = ""
        if "NegIonizable" in property_list.keys():
            for l in tuple(set(property_list["NegIonizable"])): line = line + l +" "
            sel_a = "((resname ARG LYS ) and (name NH* NZ)) or ((resname HI2 ) and (name HD HE))"
            sel_b = '((resname '+sel_ligands+") and (name "+line+") )"
            IFP_prop_list.append(IFP_prop("IN",line,sel_a,sel_b,4.5))
    except:
        print("IN failed")
        pass
    
    try:  #--- cation-aril interactions  and aromatic stacking
        line = ""
        if "PosIonizable" in property_list.keys():
            for l in tuple(set(property_list["PosIonizable"])): line = line + l +" "
        if "Aromatic" in property_list.keys():  # pi-pi
            for l in np.asarray(property_list["Aromatic"]): line = line + l +" "
        sel_b = '((resname '+sel_ligands+" ) and (name "+line+"))"
        sel_a = "((resname PHE TRP TYR HIS HIE HID HE2) and (name CZ* CD* CE* CG* CH* NE* ND*))"      
        if("PosIonizable" in property_list.keys()): 
            IFP_prop_list.append(IFP_prop("AR",line,sel_a,sel_b,4.5))
        if("Aromatic" in property_list.keys()):
            IFP_prop_list.append(IFP_prop("AR",line,sel_a,sel_b,5.5))
    except:
        print("AR1 failed")
        pass

    try: #--- S -aromatic 
        sel_b = '((resname '+sel_ligands+") and (type S ))"
        sel_a = "((resname PHE TRP TYR HI2 HIS HIE HID) and (name CZ* CD* CE* CG* CH* NE* ND*)) "
        IFP_prop_list.append(IFP_prop("AR",line,sel_a,sel_b,4.5))
    except:
        print("AR2 failed")
        pass
    
    try: #--- aromatic ligand - cation , amide, or S interactions
        line = ""
        if "Aromatic" in property_list.keys():
            for l in np.asarray(property_list["Aromatic"]): line = line + l +" "
            sel_b = '((resname '+sel_ligands+" ) and (name "+line+") )"
            sel_a = "((resname ARG LYS ) and (name NH* NZ*)) or (backbone and name H)" #  or ((type S) and (resname MET CYS))
            IFP_prop_list.append(IFP_prop("AR",line,sel_a,sel_b,4.5))
    except:
        print("AR3 failed")
        pass
    
    try: #--- halogen bonds with atromatic or backbone carbonyl oxygen
        sel_b = '((resname '+sel_ligands+" ) and ( type I CL BR Br Cl) )"
        sel_a = "((resname PHE TRP TYR HIS HIE HID ) and (name CZ* CD* CE* CG* CH* NE* ND*)) or (backbone and name O) or ((resname ASP GLU) and (name OE* OD*))  or ((resname CYS MET) and (type S))"
        IFP_prop_list.append(IFP_prop("HL","HL",sel_a,sel_b,3.5))
    except:
        print("HL failed")
        pass
    try: # water shell  
#        sel_a = "((sphzone 12.0 resname "+sel_ligands+") and (resname WAT and type O))"
        sel_a = " (resname WAT HOH SOL and type O) "
        sel_b = '(resname '+sel_ligands+" ) and ( not type H )"
        IFP_prop_list.append(IFP_prop("WA","HE",sel_a,sel_b,3.5))
    except:
        print("WA failed")
        pass
    
    if RE:
        try: # any protein-ligand contacts 
            sel_a = "(protein) and (not type H)"
            sel_b = '(resname '+sel_ligands+" ) and ( not type H )"
            IFP_prop_list.append(IFP_prop("RE","HE",sel_a,sel_b,3.5))
        except:
            print("RE failed")
            pass
        
    if len(Lipids) > 0:
        try: # any lipid-ligand contacts 
            line = ""
            for l in Lipids: line = line + l +" "
            sel_a = '((resname '+line+' ) and (not type H)) '
            sel_b = '(resname '+sel_ligands+") and (not type H)"
            IFP_prop_list.append(IFP_prop("LL",line,sel_a,sel_b,4.0))
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
def IFP(u_mem,sel_ligands,property_list, WB_analysis = True, RE = True,Lipids = []):
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
    
    h = hb.HydrogenBondAnalysis(u_mem, selection1 ='resname '+sel_ligands,selection2=' not resname WAT HOH SOL '+sel_ligands, distance=3.3, angle=110, forcefield='OtherFF')
    print("Start HB analysis",datetime.datetime.now().time())
    h.run()
    h.generate_table()
    df_HB = pd.DataFrame.from_records(h.table)
    
    #---------------------------------------------------------------
    #------ find water bridges between ligand and protein--------
    df_WB = pd.DataFrame()
    if WB_analysis :
        print("Start WB analysis",datetime.datetime.now().time())
        try:
            # will add O atom of water as an HB donor (it is concidered as acceptor by default assuming as a protein backbone atom)
            hb.WaterBridgeAnalysis.DEFAULT_DONORS['OtherFF'] = hb.WaterBridgeAnalysis.DEFAULT_DONORS['OtherFF']+tuple(set("O"))        
#            hb.WaterBridgeAnalysis.DEFAULT_ACCEPTORS['OtherFF'] = hb.WaterBridgeAnalysis.DEFAULT_DONORS['OtherFF']+tuple(set("O"))        
#            w = hb.WaterBridgeAnalysis(u_mem, 'resname '+sel_ligands, ' not resname WAT HOH SOL ',water_selection=" resname WAT HOH SOL ", 
#                distance=3.5, angle=110, forcefield='OtherFF',output_format="donor_acceptor",order=3)
            w = hb.WaterBridgeAnalysis(u_mem, selection1  = 'resname '+sel_ligands, selection2  = ' not resname WAT HOH SOL '+sel_ligands,water_selection=" resname WAT HOH SOL ",  distance=3.5, angle=110, forcefield='OtherFF',order=5)
            w.run()
            w.generate_table()
            df_WB = pd.DataFrame.from_records(w.table)
        except:
            print("Water bridge analysis failed ")
        
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
                        if(u.resid in ar_resid[ar_n > 4]): 
                            # check also residue name to deal the case of residues with the same id
                            if( np.unique(np.asarray(u_ar_n)[np.where(u_ar==u.resid)[0]]).shape[0])== 1: 
                                found.append([IFP_type.name+"_"+u.resname+str(u.resid),u.name])
 #                               print("!!!!!!!!!!!!!  pi-pi stacking found",ar_n,u.resname,str(u.resid))
                        # here we will check if cation (LYS or ARG) really contact an aromatic ring of the ligand
                        elif(u.resname in ["LYS","ARG"]):
                            if u.resname == "LYS" : cation = "LYS"
                            else: cation = "ARG"
         #                   cation = u_resname 
                            if("Aromatic" in property_list.keys()):
                                line_ar = ""
                                for l in np.asarray(property_list["Aromatic"]): line_ar = line_ar + l +" "
                                line1 = "(resname "+sel_ligands+" and ( not type H O) and name "+line_ar+") and around 4.5 (resid "+str(u.resid[u.resname == cation][0]) + " and type N)" 
                                u1_list = (u_mem.select_atoms(line1,updating=True))
                                if(len(u1_list) > 4): 
                                    found.append([IFP_type.name+"_"+u.resname+str(u.resid),u.name])
                 #               print("!!!!!!!!!!!!!  Cat-Ar  interactions found",len(u1_list),u1_list)
                        # now we check if aromatc residue (HIS) is perpendicular to the aromatic fragment of the ligand
                        ### TOBE checked if this works!!!!!================================================
                        elif(u.resname in ["PHE", "TRP", "TYR","HIS","HIE","HID","HI2"]) and (u.resid in ar_resid[ar_n < 4]):
                            if("Aromatic" in property_list.keys()):
                                line_ar = ""
                                for l in np.asarray(property_list["Aromatic"]): line_ar = line_ar + l +" "
                                line1 = "(resname "+sel_ligands+" and ( not type H O) and name "+line_ar+") and around 5. (resid "+str(u.resid) + " and (name NE* ND* CE* CD* CZ* CH*))" 
                                u1_list = (u_mem.select_atoms(line1,updating=True))
                                if(len(u1_list) > 4): 
                                    found.append([IFP_type.name+"_"+u.resname+str(u.resid),u.name])
                        ### TOBE checked if this works!!!!!================================================
                                
            else:  
                for u in u_list:
                        found.append([IFP_type.name+"_"+u.resname+str(u.resid),u.name]) 
                        if IFP_type.name == "HL": print("HAL:",u.resname,u.resid,u.name)
                
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
        df_WB_INH = df_WB[(df_WB.sele1_resnm.isin([ligand_name]) & df_WB.sele2_resnm.isin(["WAT","HOH","SOL"]))]
#        print(df_WB_INH)
        # get a list of WAT-Prot (sele1 - sele2)
        df_WB_Prot = df_WB[(~(df_WB.sele2_resnm.isin([ligand_name,"WAT","HOH","SOL"])) & (df_WB.sele1_resnm.isin(["WAT","HOH","SOL"])))]
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
def Plot_IFP(df,contact_collection=None,out_name=""):
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
    ifp_list = ["HY","AR","HD","HA","HL","IP","IN","WB","LIP"]
    columns_IFP = []  # standard IFP
    columns_CONT = []  # just contacts
    for c in df.columns.tolist():
        if c[0:2] in ifp_list:
            columns_IFP.append(c)
        elif c[0:2]  == "RE":
            columns_CONT.append(c)
                
    fig = plt.figure(figsize=(16, 6))
    gs = gridspec.GridSpec(1, 3, width_ratios=[4, 2, 1]) 
    ax = plt.subplot(gs[0])
   
    ax =  plt.subplot(gs[0])
    ax.set_title('IFP')
    df1 = df[columns_IFP].values
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
    if out_name == "":   plt.show()
    else: plt.savefig(out_name,dpi=300)
    return

########################################
#
#     get ligand chemical properties
#
########################################
def  ligand_properties(ligand_pdb,ligand_mol2):
    """
    Parameters:
    ligand_pdb - ligand structure file  in the PDB format
    ligand_mol2 - ligand structure file  in the Mol2 format (not all mol2 format work, but generated by MOE does)
    
    Returns:
    """
    ff=open(ligand_pdb,"r")
    lines = ff.readlines()
    ff.close()
    list_labels = []
    for line in lines:
        if (line.split()[0] == 'ATOM' or line.split()[0] == 'HETATM'): list_labels.append(line.split()[2]) 

    mol = Chem.rdmolfiles.MolFromMol2File(ligand_mol2,removeHs=False)   
    fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
    factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
    feats = factory.GetFeaturesForMol(mol)

    properties_list = {}
    for f in feats:
        prop = f.GetFamily()  #  get property name
        at_indx  = list(f.GetAtomIds())  # get atom index
        if prop not in properties_list.keys():
            properties_list[prop]=[]
        if(len(at_indx) > 0 ):
            for l in at_indx:
                try:
                    properties_list[prop].append(list_labels[l])
                except:
                    print("Error:",prop," atom index: ",at_indx, " are not found") 
        else: properties_list[prop].append(list_labels[at_indx[0]])
    return(properties_list)





