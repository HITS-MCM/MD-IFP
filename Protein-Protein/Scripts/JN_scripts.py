#!/usr/bin/env python
# coding: utf-8

### Collection of scripts for analysing pre-computed PP-REs dataframes from RAMD trajectories
#   to compute residence time and dissociation routes
# 
# 
#############################
### v 1.1
#    Copyright (c) 2024
#    Released under the EUPL Licence, v1.2 or any higher version
#    
### Authors: Daria Kokh, Giulia D'Arrigo
#    Heidelberg Institute of Theoretical Studies (HITS, www.h-its.org)
#    Heidelberg, Germany
#
### Contact: mcmsoft.h-its.org 
# 
### Packages required:
#     numpy
#     pandas
#     matplotlib
#     seaborn
#     MDAnalysis
#     code is written on Python 3 and tested on the version 3.7
##########################################################################


import glob, os
import sys
import subprocess
import numpy as np

import pandas as pd
from pandas import ExcelFile 

from matplotlib import *
import matplotlib as mpl
from matplotlib import cm
from matplotlib.patches import Arrow, Circle
import matplotlib.gridspec as GS
import matplotlib.ticker
import  pylab as plt
import seaborn
import seaborn as sns


#import ipywidgets as widgets
from sklearn import linear_model
from sklearn import preprocessing
from sklearn.cluster import KMeans
from matplotlib import gridspec
from sklearn.metrics import r2_score

from scipy import stats
from scipy.optimize import curve_fit
from scipy.stats import norm

##############################################
#
#Bootstrapping analysis
#
##############################################

def bootstrapp(t, rounds=50000):
    max_shuffle = rounds
    alpha = 0.8
    sub_set = int(alpha*len(t))
    tau_bootstr = []
    for i in range(1,max_shuffle):
        # generate a sub-set
        np.random.shuffle(t)
        t_b = t[:sub_set]
        # find residence time from a sub-stet
        t_b_sorted_50 =(np.sort(t_b)[int(len(t_b)/2.0-0.5)]+np.sort(t_b)[int(len(t_b)/2)])/2.0
        tau_bootstr.append(t_b_sorted_50)
    return(tau_bootstr)

##############################################
#
#Different criteria of residence time
#
##############################################

def fun_few_contact_first(k,cols,dock_resi,threshold,r,tr):
    k1 = k[cols]
    if k1[k1 < threshold*10].shape[0] > 0:
        contacts = k1[k1 < threshold*10].count(axis = 1)
        if k.time.values[(contacts < int(len(dock_resi)/2) ).values].shape[0] > 0:
            time = k.time.values[(contacts < int(len(dock_resi)/2) ).values][0]
        else:
            print("Traj ",r,tr,"Few_contacts_first not found: no dissociation? ") 
            time = k.time.values[-1]
    else:
        print("Traj ",r,tr,"Few_contacts_first not found: no dissociation? " )
        time = k.time.values[-1]
    return (time)

def fun_many_contact_last(k,cols,dock_resi,threshold,r,tr):
    k1 = k[cols]
    contacts = k1[k1 < threshold*10].count(axis = 1)
    if k.time.values[(contacts > int(len(dock_resi)/2))].shape[0] > 0:
        time = k.time.values[(contacts > int(len(dock_resi)/2)).values][-1]
    else:
        print(r,tr,"Many_contacts_last not found: no dissociation ")
        time =k.time.values[1]
    return (time)

def fun_resi_last(k,dock_resi,threshold,r,tr):
    k1 = k[dock_resi]
    k2 = k.time.values[k1.values.mean(axis=1) < threshold*10 ]
    if len(k2) > 1:
        time =k2[-1]
    else:
        print("Traj ",r,tr,"By_residue_last not found: no proper binding-site contacts? " )
        time =k.time.values[1]
    return(time)

def fun_resi_first(k,dock_resi,threshold,r,tr):
    k1 = k[dock_resi]  
    k2 = k.time.values[k1.values.mean(axis=1) > threshold*10 ] 
    if len(k2) > 1:
        time = k2[1]
    else:
        print("Traj ",r,tr,"The first unbinding event was not found: no dissociation? " )
        time = k.time.values[-1]
    return(time)

def fun_resi_last1(k,dock_resi,threshold,r,tr):
    k1 = k[dock_resi]
    contacts = k1[k1 < threshold*10].count(axis = 1)
    try:
        time = k.time.values[(contacts > 5).values][-1]
    except:
        print(r,tr,"By_residue_last not found: no proper binding-site contacts? " )
        time =k.time.values[1]
    return (time)

def fun_resi_first1(k,dock_resi,threshold,r,tr):
    k1 = k[dock_resi]  
    contacts = k1[k1 < threshold*10].count(axis = 1)
    try:
        contacts = k1[k1 < threshold*10].count(axis = 1)
        time = k.time.values[(contacts < 5).values][0]
    except:
        print("Traj ",r,tr,"The first unbinding event was not found: no dissociation? " )
        time = k.time.values[-1]
    return (time)

def fun_nocontact(k,dock_resi,threshold,r,tr):
    k1 = k[dock_resi]  
    contacts = k1[k1 < threshold*10].count(axis = 1)
    try:
        contacts = k1[k1 < threshold*10].count(axis = 1)
        time = k.time.values[(contacts ==0).values][0]
    except:
        print("Traj ",r,tr,"contact 0 not found: no dissociation? " )
        time = k.time.values[-1]
    return (time)

##############################################
#
#The following two functions are to plot the 
#distribution of residence time for each criterion
#
##############################################

def tau_diss_event(df_wt,cols,dock_resi,threshold =5.5, tau_list_name = ["COM-COM","By_residue_last","Many_contacts_last","Few_contacts_first","By_residue_first"]):
    """
    Parameters:
    Results:
    """
    
    tau_list = {}
    tau_list_replT = []
    tau_list_wT = []
    tau_list_cT = []
    tau_list_rT = []
    tau_list_cT1 = []
    tau_list_rT1 = []
    tau_list_nT1 =[]
    tau_list_rT2 = []
    print(np.unique(df_wt.Repl))
    for r in np.unique(df_wt.Repl):
        tau = []
        df_r = df_wt[df_wt.Repl == r]
        if "Water" in tau_list_name : 
            dfW_r = dfW_wt[dfW_wt.Repl == r]
        tau_list_repl = []
        tau_list_w = []
        tau_list_c = []
        tau_list_r = []
        tau_list_c1 = []
        tau_list_r1 = []
        tau_list_n1 =[]
        tau_list_r2 = []
        for ti, tr in enumerate(np.unique(df_r.Traj)):
            k = df_r[df_r.Traj == tr]
            if k.time.values.shape[0] < 2: continue
            # standard COM-COM time - by threshold o COM-COM distance
            if "COM-COM" in tau_list_name : 
                tau_list_repl.append(k.time.values[-1])
            # by the number of interfacial water (the first snapshot where there is only 10 molecules of  interfacial water - dissociation)
            if "Water" in tau_list_name : 
                kW = dfW_r[dfW_r.Traj == tr]
                tau_list_w.append(kW[kW.Water < 10].time.values[0])

            #  by number of contacts below a threshold distance (No_contacts - dissociation) 
            if "Few_contacts_first" in tau_list_name : 
                tau_list_c.append(fun_few_contact_first(k,cols,dock_resi,threshold,r,tr))
            if len(dock_resi) > 0:
                # by the average distance between binding site residues ( the last snapshot below a threshold - dissociation)
                if "By_residue_last" in tau_list_name :
                    tau_list_r.append(fun_resi_last(k,dock_resi,threshold,r,tr))
                # by the average distance between binding site residues ( the first snapshot below a threshold - dissociation)
                if "By_residue_first" in tau_list_name :
                    tau_list_r1.append(fun_resi_first(k,dock_resi,threshold,r,tr))
                if "By_residue_first1" in tau_list_name :
                    tau_list_r2.append(fun_resi_first1(k,dock_resi,threshold,r,tr))
                # by the number of binding site residues with the distance above a threshold ( the first snapshot with more than 50% - dissociation)
                if "Many_contacts_last" in tau_list_name : 
                    tau_list_c1.append(fun_many_contact_last(k,cols,dock_resi,threshold,r,tr))
                if "No_contact" in tau_list_name :
                    tau_list_n1 =(fun_nocontact(k,dock_resi,threshold,r,tr))
                    
                    
        if "COM-COM" in tau_list_name :         tau_list_replT.append(tau_list_repl)
        if "Water" in tau_list_name :          tau_list_wT.append(tau_list_w)
        if "Few_contacts_first" in tau_list_name :       tau_list_cT.append(tau_list_c)
        if "By_residue_last" in tau_list_name :         tau_list_rT.append(tau_list_r)
        if "Many_contacts_last" in tau_list_name :       tau_list_cT1.append(tau_list_c1)
        if "By_residue_first" in tau_list_name :         tau_list_rT1.append(tau_list_r1)
        if "No_contact" in tau_list_name :         tau_list_nT1.append(tau_list_n1)
        if "By_residue_first1" in tau_list_name :         tau_list_rT2.append(tau_list_r2)

    if "COM-COM" in tau_list_name :      tau_list["COM-COM"]=tau_list_replT
    if "Water" in tau_list_name :       tau_list["Water"]=tau_list_wT
    if "Few_contacts_first" in tau_list_name :     tau_list["Few_contacts_first"]=tau_list_cT
    if "By_residue_last" in tau_list_name :       tau_list["By_residue_last"]=tau_list_rT
    if "Many_contacts_last" in tau_list_name :     tau_list["Many_contacts_last"]=tau_list_cT1
    if "By_residue_first" in tau_list_name :       tau_list["By_residue_first"]=tau_list_rT1
    if "No_contact" in tau_list_name :       tau_list["No_contact"]=tau_list_nT1
    if "By_residue_first1" in tau_list_name :         tau_list_rT1.append(tau_list_r2)   
        
    #print(tau_list)
    return (tau_list)

def plot_hist_tau_list(df_D_mu,name_mu,bond_contacts,cols_mu,lim,threshold, tau_list_name = ["COM-COM","Few_contacts_first","By_residue_last","No_contact"]):
    """
    Parameters:
    Results:
    """
    bins = 10
    meanpointprops = dict(linestyle='--', linewidth=1.5, color='firebrick')
    medianpointprops = dict(linestyle='-', linewidth=2.0, color='orange')
    # loop over methods
    tau_mu = []
    #loop over mutants
    for i,(df,name,bc,cols) in enumerate(zip(df_D_mu,name_mu,bond_contacts,cols_mu)): 
        print("------------",name)
        repls= df.Repl.unique()
        tau_list=tau_diss_event(df,cols,bc,threshold =threshold,tau_list_name = tau_list_name)
        fig = plt.figure(figsize = (3*len(tau_list.values()),1),facecolor='w',dpi=150)
        gs = gridspec.GridSpec(1,len(tau_list_name), wspace=0.5) 
        ik = 0
        tau_results = {}
        # loop over methods

        for t_r,name_method in zip(tau_list.values(),tau_list.keys()):
            # print(name_method,len(t_r),len(t_r[0]))
            tau_r = []
            all_t = []
            ax1 = plt.subplot(gs[ik])
            # loop over replicas
            for (rr,t) in zip(repls,t_r):
                if len(t) >= 8:
                    for k in range(len(t),15): t.append(40*100)
                    tau = bootstrapp(np.array(t)/100., rounds=5000)
                    all_t.append(tau)
                    h, bins = np.histogram(tau,range=lim,bins=bins, density = True)
 #                   mu, std = norm.fit(np.array(t)/100.)
                    tau_r.append(np.mean(tau))
                    print (name_method,rr, t)
                else:
                    print(name_method,"replica: ",rr,": only "+str(len(t))+" trajectories")
                    continue
          #      print("---------------- ",name_method,rr,np.array(t)/100,np.mean(tau),np.std(tau))
            tau_results[name_method]=(np.mean(np.array(tau_r)),np.std(np.array(tau_r)))
            plt.title(name_method+" "+str(np.round(np.mean(np.array(tau_r)),2))+"+-"+str(np.round(np.std(np.array(tau_r)),2)),fontsize=8)
            ax1.boxplot(all_t,showmeans=True, meanline=True,meanprops=meanpointprops,medianprops = medianpointprops, bootstrap=5000) #labels = mue_set)
            #ax1.set_xticklabels(repls)
            plt.grid(color='gray', linestyle='--', alpha=0.5,linewidth=0.5)
            ax1.set_ylim(lim[0],lim[1])
            ax1.set_xlim(1,10)
            #plt.xticks(range(1,10), fontsize=8)
            plt.yticks(range(lim[0],lim[1],5), fontsize=8)
            ik += 1
        tau_mu.append(tau_results)
        plt.show()
    return(tau_mu)

##############################################
#
#Function to calculate the errors for the computed and experimental values 
#
##############################################

def Residence_time (method_names, exp_data, name_mu, tau_mu):
    
    
    method_tau = []
    method_tau_SD = []
    method_tau_err_p = []
    method_tau_err_n = []
    name_prot = []
    
    for i,method in enumerate(method_names): 
        method_tau.append([])
        method_tau_SD.append([])
        method_tau_err_p.append([])
        method_tau_err_n.append([])

    exp_list = []
    exp_list_std = []
    exp_list1 = []
    
    for i,(name,mut) in enumerate(zip(name_mu, tau_mu)):
        if (name in exp_data.prot.values) or (name[:-1] in exp_data.prot.values):
            exp_list.append(exp_data[exp_data.prot.isin([name,name[:-1]])]["koff"].values[0]) 
#             print(exp_list)
            exp_list_std.append(exp_data[exp_data.prot.isin([name,name[:-1]])]["SD_koff"].values[0])
#             print(exp_list_std)
            for j,method in enumerate(method_names): 
                method_tau[j].append(mut[method][0])
                method_tau_SD[j].append(mut[method][1])
        else:
            print("exp data for the mut "+name+" are missing")
        name_prot.append(name)
       
    for j,method in enumerate(method_names):
        method_tau_err_p[j] = (np.log10(np.array(method_tau[j])+np.array(method_tau_SD[j])))- np.log10(np.array(method_tau[j]))
        if np.isnan(method_tau_err_p[j][-1]):
            method_tau_err_p[j][-1]=0.0
        method_tau_err_n[j] = np.log10(np.array(method_tau[j])) -(np.log10(np.array(method_tau[j])-np.array(method_tau_SD[j])))
        if np.isnan(method_tau_err_n[j][-1]):
            print("yes")
            method_tau_err_n[j][-1]=0.0
#         print(method_names[j], method_tau_err_n[j])
        
        
    exp_list_std = np.array(exp_list_std)
    #print(exp_list_std)
    exp_list = np.array(exp_list)
#     print(exp_list)
    exp = np.log10(1/exp_list)
    exp_err_n = -np.log10(1/(exp_list+exp_list_std))+ exp
    exp_err_p = (np.log10(1/(exp_list-exp_list_std)))-exp
    return(exp,exp_err_n,exp_err_p,name_prot,method_tau,method_tau_err_n,method_tau_err_p)

##############################################
#
#Function to calculate the errors for the computed and experimental values 
#
##############################################

def Residence_time_kd(method_names, exp_data, name_mu, tau_mu):
    
    
    method_tau = []
    method_tau_SD = []
    method_tau_err_p = []
    method_tau_err_n = []
    name_prot = []
    
    for i,method in enumerate(method_names): 
        method_tau.append([])
        method_tau_SD.append([])
        method_tau_err_p.append([])
        method_tau_err_n.append([])

    exp_list_kd = []
    exp_list_std_kd = []
    exp_list1 = []
    
    for i,(name,mut) in enumerate(zip(name_mu, tau_mu)):
        if (name in exp_data.prot.values) or (name[:-1] in exp_data.prot.values):
            exp_list_kd.append(exp_data[exp_data.prot.isin([name,name[:-1]])]["KD"].values[0]) 
#             print(exp_list)
            exp_list_std_kd.append(exp_data[exp_data.prot.isin([name,name[:-1]])]["SD_KD"].values[0])
#             print(exp_list_std)
            for j,method in enumerate(method_names): 
                method_tau[j].append(mut[method][0])
                method_tau_SD[j].append(mut[method][1])
        else:
            print("exp data for the mut "+name+" are missing")
        name_prot.append(name)
       
    for j,method in enumerate(method_names):
        method_tau_err_p[j] = (np.log10(np.array(method_tau[j])+np.array(method_tau_SD[j])))- np.log10(np.array(method_tau[j]))
        if np.isnan(method_tau_err_p[j][-1]):
            method_tau_err_p[j][-1]=0.0
        method_tau_err_n[j] = np.log10(np.array(method_tau[j])) -(np.log10(np.array(method_tau[j])-np.array(method_tau_SD[j])))
        if np.isnan(method_tau_err_n[j][-1]):
            print("yes")
            method_tau_err_n[j][-1]=0.0
#         print(method_names[j], method_tau_err_n[j])
        
        
    exp_list_std_kd = np.array(exp_list_std_kd)
#     print(exp_liststd)
    exp_list_kd = np.array(exp_list_kd)
    #print(exp_list_kd)
    exp_kd = np.log10(1/exp_list_kd)
    #print(exp_kd)
    exp_err_n_kd = (np.log10(1/(exp_list_kd-exp_list_std_kd)))-exp_kd
    exp_err_p_kd = (np.log10(1/(exp_list_kd-exp_list_std_kd)))-exp_kd
    return(exp_kd,exp_err_n_kd,exp_err_p_kd,name_prot,method_tau,method_tau_err_n,method_tau_err_p)

##############################################
#
#Function to compute RMSD of equilibration trajectories
#
##############################################

def plot_RMSD_eq(df_E_mu,name_mu):

    boxprops = dict(linestyle='-', linewidth=0.5, color='k')
    whiskerprops= dict(linestyle='-', linewidth=0.5, color='k')
    capprops=dict(linestyle='-', linewidth=0.5, color='k')
    flierprops = dict(markerfacecolor='k', marker='o', markersize=1)
    medianprops = dict(linestyle='-', linewidth=1, color='black')
    meanpointprops = dict(marker='o',markersize=1.5, markeredgecolor='purple',
                      markerfacecolor='purple')

    # loop over methods
    ik = 0
    #loop over mutants
    fig = plt.figure(figsize = (8,10),facecolor='w',dpi=300)
    gs = gridspec.GridSpec(nrows=6, ncols=4, wspace=0.5,hspace=0.7)
    for i,(df,name) in enumerate(zip(df_E_mu,name_mu)):
        repls = df.Repl.unique()
        all_t = []
        for j,r in enumerate(repls):
            all_t.append(df[df.Repl==r].RMSD12.values)
        ax1 = plt.subplot(gs[ik])
        plt.title(name,fontsize=8)
        ax1.boxplot(all_t,showmeans=True, meanline=False,patch_artist=False,notch=False, 
            bootstrap=5000,boxprops=boxprops,medianprops=medianprops,whiskerprops=whiskerprops,
            flierprops=flierprops, capprops=capprops,meanprops=meanpointprops)
        ax1.set_xticklabels(repls,fontsize=7.5)
        ax1.set_yticks(range(1,6))
        ax1.set_yticklabels(range(1,6),fontsize=7.5)
        plt.grid(color='gray', linestyle='--', alpha=0.3,linewidth=0.3)
        ax1.set_ylim(0.5,6)
        ax1.set_xlim(0.5,len(repls)+0.5)
        plt.ylabel("RMSD (Å)",fontsize=8)
        ik += 1
    return

##############################################
#
#Function to compute interfacial waters
#
##############################################

def mean_repl_Water_inter(dfW_wt_list,name_list):
    
    
    boxprops = dict(linestyle='-', linewidth=0.5, color='k')
    whiskerprops= dict(linestyle='-', linewidth=0.5, color='k')
    capprops=dict(linestyle='-', linewidth=0.5, color='k')
    flierprops = dict(markerfacecolor='k', marker='o', markersize=1)
    medianprops = dict(linestyle='-', linewidth=1, color='black')
    meanpointprops = dict(marker='o',markersize=1.5, markeredgecolor='purple',
                      markerfacecolor='purple')

    fig = plt.figure(figsize = (8,10),facecolor='w',dpi=300)
    gs = gridspec.GridSpec(nrows=6, ncols=4, wspace=0.5,hspace=0.7)

    for i,(dfW_wt,name) in enumerate(zip(dfW_wt_list,name_list)):
        w = []
        wt = []
        ax1 = plt.subplot(gs[i])
#         max_wat = 10*(int(np.max(dfW_wt.Water.values)/10+1))
        repls = np.unique(dfW_wt.Repl.values).astype(int)
        ax1.set_ylim(10,45)
        ax1.set_xlim(0.5,len(repls)+0.5)
        ax1.grid(color='gray', linestyle='--', alpha=0.3,linewidth=0.3)
        plt.ylabel("Interfacial Waters",fontsize=8)
        plt.title(name,fontsize=8)
#         ax1.set_yticks(range(0, 45))
        ax1.yaxis.set_ticks(np.arange(10, 45, 10))
        ax1.set_xticklabels(repls,fontsize=7.5)
        ax1.set_yticklabels(np.arange(10,45,10),fontsize=7.5)
        for t in np.unique(dfW_wt.Traj.values):
            w.append(dfW_wt[dfW_wt.Traj == t].Water.values)
            wt.append(dfW_wt[dfW_wt.Traj == t].WaterTot.values)
        ax1.boxplot(wt,showmeans=True, meanline=False,patch_artist=False,notch=False, 
            bootstrap=5000,boxprops=boxprops,medianprops=medianprops,whiskerprops=whiskerprops,
            flierprops=flierprops, capprops=capprops,meanprops=meanpointprops)
#         ax1.violinplot(w, showmeans=False,showmedians=True, showextrema=True,widths=0.5)
#     plt.xlabel("Replica",fontsize=7)
    # fig.show()
    return 

##############################################
#
#Function to compute buried waters
#
##############################################

def mean_repl_Water_buried(dfW_wt_list,name_list):
    
    
    boxprops = dict(linestyle='-', linewidth=0.5, color='k')
    whiskerprops= dict(linestyle='-', linewidth=0.5, color='k')
    capprops=dict(linestyle='-', linewidth=0.5, color='k')
    flierprops = dict(markerfacecolor='k', marker='o', markersize=1)
    medianprops = dict(linestyle='-', linewidth=1, color='black')
    meanpointprops = dict(marker='o',markersize=1.5, markeredgecolor='purple',
                      markerfacecolor='purple')

    fig = plt.figure(figsize = (8,10),facecolor='w',dpi=300)
    gs = gridspec.GridSpec(nrows=6, ncols=4, wspace=0.5,hspace=0.7)

    for i,(dfW_wt,name) in enumerate(zip(dfW_wt_list,name_list)):
        w = []
        wt = []
        ax1 = plt.subplot(gs[i])
#         max_wat = 10*(int(np.max(dfW_wt.Water.values)/10+1))
        repls = np.unique(dfW_wt.Repl.values).astype(int)
        ax1.set_ylim(5,30)
        ax1.set_xlim(0.5,len(repls)+0.5)
        ax1.grid(color='gray', linestyle='--', alpha=0.3,linewidth=0.3)
        plt.ylabel("Buried Waters",fontsize=8)
        plt.title(name,fontsize=8)
#         ax1.set_yticks(range(0, 45))
        ax1.yaxis.set_ticks(np.arange(5, 30, 10))
        ax1.set_xticklabels(repls,fontsize=7.5)
        ax1.set_yticklabels(np.arange(5,30,10),fontsize=7.5)
        for t in np.unique(dfW_wt.Traj.values):
            w.append(dfW_wt[dfW_wt.Traj == t].Water.values)
            wt.append(dfW_wt[dfW_wt.Traj == t].WaterTot.values)
        ax1.boxplot(w,showmeans=True, meanline=False,patch_artist=False,notch=False, 
            bootstrap=5000,boxprops=boxprops,medianprops=medianprops,whiskerprops=whiskerprops,
            flierprops=flierprops, capprops=capprops,meanprops=meanpointprops)
    return

##############################################
#
#Graph-based representation of protein-protein dissociation trajectories
#
##############################################

def plot_graph_COM(df_ext,file_save = "",ligand = "",draw_round = True,water = False,edges_show=True,diss_show = True):

    """  
    Parameters:
    df_ext - IFP database
    file_save - file name to save an image
    ligand - generate dissociation pathways for a selected ligand only (note, that clusters properties still include data for all ligands)
    draw_round - type of representation (plain/round)
    water - visualize number of water molecules in the ligand solvation shell for each clutser
    
    Returns:
    cluster label sorted by increase of the average  RMSD in each cluster
    """

    from matplotlib.patches import ArrowStyle
    from matplotlib.patches import Ellipse
    from scipy.spatial import distance
    df_ext_ligand = df_ext
    if len(ligand)> 0:
        try:
            df_ext_ligand = df_ext[df_ext.ligand.isin(ligand)]
            print("Edges and node size will be shown for one ligand:",ligand)
        except:
            print("ligand "+ligand+" was not found in the database. Whole database will be analyzed")
   
    #-------- collect data ------------------------
    label_rmsd = []  
    label_com = []   
    label_size = []
    label_water = []

    labels_list,nodes = np.unique(df_ext.label.values,return_counts= True)
    edges = np.zeros((labels_list.shape[0],labels_list.shape[0]),dtype = float)
    coms = np.zeros((labels_list.shape[0],labels_list.shape[0]),dtype = float)
    label_egress = np.zeros((labels_list.shape[0]),dtype = float)

    df_min_rmsd = df_ext[df_ext.RMSDl == df_ext.RMSDl.min()]
    min_rmsd = df_min_rmsd.RMSDl.values[0]
    min_COM = df_min_rmsd.COM.values[0]
    
    # cluster properies: firxt loop over clusters
    for i,l in enumerate(labels_list):
        t_lig = df_ext_ligand[df_ext_ligand.label == l]
        t = df_ext[df_ext.label == l]
        label_rmsd.append(t.RMSDl.mean())
        
        com_av = []
        for com in t.COM.values:
            com_av.append(max(0,com-min_COM))
        label_com.append(np.mean(com_av))
        if t_lig.shape[0] > 0:
#             print("number frames in each cluster", l,t_lig.shape[0], "number total frames",df_ext_ligand.shape[0]) #how many frames in each cluster
            label_size.append(100*t.shape[0]/df_ext.shape[0]) 
#             print(label_size)
        else:
            label_size.append(0)
        if water:
            label_water.append(int(t.WaterTot.mean()))
        for j in range(0,i):
            coms[i,j] = distance.euclidean(label_com[i],label_com[j])
#         print("cluster",l, "size",label_size,t.shape[0],df_ext_ligand.shape[0])
#         print("cluster ",l,"STD of COM: ", t.COM.std(),"STD of RMSD: ",t.RMSDl.std(),"Water:",label_water)
    
    # transitions between clusters:
    for l,(df_label,df_time) in enumerate(zip(df_ext_ligand.label.values,df_ext_ligand.time.values)):
        thiselem = df_ext_ligand.Traj.values[l-1]
        nextelem = df_ext_ligand.Traj.values[(l) % len(df_ext_ligand.Traj.values)]
        if df_time != 0   and l != 0:
            if(thiselem == nextelem) and (df_ext_ligand.label.values[l-1] != df_label):
                edges[df_ext_ligand.label.values[l-1],df_label] += labels_list.shape[0]/df_ext_ligand.label.values.shape[0]
                edges_df= pd.DataFrame(edges, columns = labels_list)
            else:
                continue  

    # find last frames in  trajectories
    for j in np.unique(df_ext.Repl.values):
        df_ext_Repl = df_ext[df_ext.Repl == j]
        for i in np.unique(df_ext_Repl.Traj.values.astype(int)):
            df_ext_Repl_Traj = df_ext_Repl[df_ext_Repl.Traj == str(i)]
            label_egress[df_ext_Repl_Traj.label.values[-1]] += 100./df_ext_ligand.label.values.shape[0] 
       
    #-----------------------------------------------------
    
    #-------- Coloring of nodes by the average RMSD in the clusters-----------
    #indx_first_com = np.argwhere((np.asarray(label_rmsd) == min(label_rmsd)))[0][0]
    dist_rmsd = []
    for i,l in enumerate(labels_list): 
        dist_rmsd.append(np.round(np.abs(label_rmsd[i]-min_rmsd),2))
    dist_com = (10*np.asarray(label_com)/np.max(label_com)).astype(int)
    dist_rmsd = (10*np.asarray(dist_rmsd)/np.max(dist_rmsd)).astype(int)
    # print("COM displacement in each cluster to be used for plotting:",dist_com)
    # print("RMSDs displacement in each cluster to be used for plotting:",dist_rmsd)
    print("COM min:",min_COM)
    print("RMSD min:",min_rmsd)
    #------------------------------------------------
    # print( min(label_rmsd),dist_com)
    fig = plt.figure(figsize = (6,2),facecolor='w',dpi=150) 
    gs = GS.GridSpec(1,2, width_ratios=[1,1],wspace=0.2) 
#     ax = plt.subplot(gs[0]) 
#     plt.title("Transition density")
#     ec = plt.imshow(edges,cmap='Oranges')
#     cbar = plt.colorbar(ec)
#     cbar.set_label('density')
#     ax = plt.subplot(gs[1])
#     plt.title("Flow")
    flow_tot = edges-edges.T
#     ec1 = plt.imshow(flow_tot,cmap='Greys')
#     cbar = plt.colorbar(ec1)
#     cbar.set_label('density')
#     plt.show()
    #print("Flow:\n",(flow*10000.).astype(int))
    
    #------------------------------------------------------------

    starting_labels = df_ext[df_ext.time == 0].label.values # list of clusters of all first frames in all trajectories
    starting_list, starting_count = np.unique(starting_labels, return_counts=True) # list of clusters that appear as first frame

    #------------ older cluster position-------------
    label2scale = label_com #label_rmsd 
    label2order = labels_list #nodes 
    index_x_t = np.argsort(label2scale)
    # x and y positions of each cluster in the plot
    label_x = np.zeros((len(labels_list)),dtype = float)  
    label_y = np.zeros((len(labels_list)),dtype = float)  
    color_com = np.zeros((len(labels_list)),dtype = float)  
    max_x = max(label2scale)
    step_x = 1.0*max_x/max(label2scale)
    # firt clusters that contain first frames of trajectories
    for i,s in enumerate(starting_list):
        label_x[s] = label2scale[s]*step_x
        label_y[s] = label2order[s]
        color_com[s] = 10*dist_rmsd[s] #dist_com[s]

    # then the rest 
    j = 0
    for l in index_x_t:
        if (labels_list[l] not in starting_list):
            while (label_x[j] != 0):
                j += 1
                if(j == len(labels_list)): break
            if(j == len(labels_list)):break
            label_x[labels_list[l]] = label2scale[l]*step_x 
            label_y[labels_list[l]] = label2order[l]
            color_com[labels_list[l]] = 10*dist_rmsd[l] #dist_com[l]

    # since the last node is usually very far from all others we will damp its color
    color_com[color_com == np.max(color_com)] = max(int(np.sort(color_com)[-2]+0.25*(np.sort(color_com)[-1]-np.sort(color_com)[-2] )),np.sort(color_com)[-1]-45) 

    # set logarythmic scale
    label_x = np.log10(label_x)
    x_tick_lable = []
    x_tick_pos = []
    for k in range(0,2):
        for ii,i in enumerate(range(pow(10,k),pow(10,k+1),pow(10,k))):  
            x_tick_lable.append(str(i))
            x_tick_pos.append(np.log10(i))
            if(i > 35): break
      
    if draw_round:
        alpha = 0.9*2*3.14*label_x/np.max(label_x)
        alpha_regular = 0.9*2*3.14*np.asarray(x_tick_pos)/max(x_tick_pos)
        label_y = np.sin(alpha)
        label_x = np.cos(alpha)
        fig = plt.figure(figsize=(6, 6),dpi=150)
        gs = GS.GridSpec(1, 1) #, width_ratios=[1, 1]) 
        ax = plt.subplot(gs[0])
        plt.scatter(x=np.cos(alpha_regular),y=np.sin(alpha_regular), c='k',s=10)
        for l,p in zip(x_tick_lable,x_tick_pos):
            ax.annotate(str(l)+"A", (1.2*np.cos(0.9*2*3.14*p/max(x_tick_pos)),np.sin(0.9*2*3.14*p/max(x_tick_pos))),fontsize=14,color="gray")
        plt.xlim(-1.3,1.3)
        plt.ylim(-1.3,1.3)
    else:
        fig = plt.figure(figsize=(3, 2),dpi=300)
        gs = GS.GridSpec(1, 1) #, width_ratios=[1, 1]) 
        ax = plt.subplot(gs[0])
#         ax.set_ylabel('Cluster', fontsize=8)
        ax.set_xlabel(r'< $\Delta{COM}  > [Angstrom]$', fontsize=7) 
        plt.xticks(x_tick_pos,x_tick_lable, fontsize=6)
        plt.yticks([])
        ax.tick_params(labelsize=6)
        plt.ylim(-0.8,len(label_y)+0.8)
        plt.xlim(-0.5,1*np.max(label_x))         
#         plt.title(ligand[0],fontsize = 6,color = "k", position=(-1, 0))
        for i,x in enumerate(np.sort(label_x)):
            plt.text(x,max(label_y)+2,str(i+1), fontsize=6)
 #       plt.grid()

    # flow and magnitude
    el = Ellipse((2, -1), 0.4, 0.4)    
    max_flow = 0
    for l in range(0,label_x.shape[0]):
        for n in range(l+1,label_x.shape[0]):
            max_flow  = max(max_flow,edges[l,n] - edges[n,l])
            
    for l in range(0,label_x.shape[0]):
        for n in range(l+1,label_x.shape[0]):
            # total number of transitions in both directions
            if (label_rmsd[l] > label_rmsd[n]):
                a = n
                b = l
            else:
                a = l
                b = n
            xy=(label_x[b],label_y[b])
            xytext=(label_x[a],label_y[a]) 
            if edges_show :
                if (edges[l,n] > 0) :
                    if  (np.abs((label_rmsd[l] - label_rmsd[n])) > 0.5* min(label_rmsd)) or (draw_round == False):
                        ax.annotate("", xy=xy, xycoords='data',
                            xytext=xytext, textcoords='data',
                            size=edges[l,n]*500,
                            arrowprops=dict(arrowstyle="Fancy,head_length=0.2, head_width=0.4, tail_width=0.2", 
                                 fc="orange", ec="none", alpha=0.4 ,
                                 connectionstyle="arc3,rad=-0.5"),
                            )
                    if  (np.abs((label_rmsd[l] - label_rmsd[n])) > 0.5* min(label_rmsd)) or (draw_round == False):
                        ax.annotate("", xy=xytext, xycoords='data',
                            xytext=xy, textcoords='data',
                            size=edges[l,n]*600,
                            arrowprops=dict(arrowstyle="Fancy,head_length=0.2, head_width=0.4, tail_width=0.2", 
                                fc="orange", ec="none", alpha=0.4 ,
                                connectionstyle="arc3,rad=-0.5"),
                            )
            #  the flow
            flow = edges[l,n] - edges[n,l]  # flow l ----> n
            if (flow > 0) :
                a = l
                b = n
            elif (flow==0):
                a=0
                b=0
            else:
                a = n
                b = l
#             print('flow',l,n,flow)
            xy=(label_x[b],label_y[b])
            xytext=(label_x[a],label_y[a])   
            ax.annotate("", xy=xy, xycoords='data',
                        xytext=xytext, textcoords='data',
                        size=np.abs(flow)*30/max_flow,
                        arrowprops=dict(arrowstyle="Fancy,head_length=0.2, head_width=0.4, tail_width=0.2", 
                                fc="0.6", ec="none", alpha=0.8 ,
                                connectionstyle="arc3,rad=-0.5"),
                        )
            
            
    # -----------------------plot nodes
    if diss_show:
        for x,y,egress in zip(label_x,label_y,label_egress/np.sum(label_egress)):
            if egress > 0:
                ax.scatter(x,y,facecolors='none',s=500*np.asarray(max(label_size)+egress/np.sum(label_egress)),color='r',lw=2)
    if water:
        ax.scatter(label_x,label_y,facecolors='none',c=color_com,edgecolors='powderblue',s=1000*np.asarray(label_size)/np.asarray(label_size)[-1],cmap='Oranges',\
               linewidths=0.3*np.asarray(label_water))
#         print("WATERS:",np.asarray(label_water))
    else:
        ax.scatter(label_x,label_y,facecolors='none',c=color_com,edgecolors="k",\
                   s=10*np.asarray(label_size),linewidth=0.7,cmap='Oranges') 
    #---------------------------------------------
    
    if file_save != "": plt.savefig(file_save,dpi=300)  
    else:    plt.show()
    return(np.argsort(label_com),flow_tot)




