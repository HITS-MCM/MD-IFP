import glob, os
import sys
import subprocess
import numpy as np

import pandas as pd
from pandas import ExcelFile 

from matplotlib import *
from matplotlib import cm
import matplotlib.ticker
import  pylab as plt
import seaborn
import seaborn as sns

#import ipywidgets as widgets

from scipy import stats



from Scripts.IFP_generation import *


from sklearn import linear_model
from sklearn import preprocessing
from sklearn.cluster import KMeans

from kmodes.kmodes import KModes


########################################
#
#  PLOT results of the clastering analysis
#  plot only trajectories where ligand completely dissociated (RMSD of the last snapshot is > 10A) 
#  plot only trajectories where ligand  started with RMSD < 10A
#
########################################

def plot_graph(df_ext,cluster_show=-1,replica_show = False):
    """
    Parameters:
    df_ext - IFP database
    cluster_show - cluster of trajectories to be highlighed 
    Returns:
    """
    """   
    label_time = []
    label_timeSD = []
    label_rmsd = []
    label_rmsdSD = []
    label_wat = []
    label_watSD = []
    label_repl = []
    label_size = []
    label_rgyr = []
    label_rgyrSD = []
    label_length = []
    label_lengthSD = []


    labels_list = np.unique(df_ext.label.tolist())
    for l in labels_list:
        t = df_ext[df_ext.label == l]
        label_time.append(int(t.time.mean()))
        if(len(t.time.tolist())>1):     label_timeSD.append(int(t.time.std()))
        else:  label_timeSD.append(0)
        label_length.append(int(t.length.mean()))
        if(len(t.time.tolist())>1):     label_lengthSD.append(int(t.length.std()))
        else:  label_lengthSD.append(0)
        label_rmsd.append(t.RMSDl.mean())
        if(len(t.time.tolist())>1):     label_rmsdSD.append(t.RMSDl.std())
        else:  label_rmsdSD.append(0)
        label_wat.append(int(t.WAT.mean()))
        if(len(t.time.tolist())>1):     label_watSD.append(int(t.WAT.std()))
        else:  label_watSD.append(0)
        label_rgyr.append(t.RGyr.mean())
        if(len(t.time.tolist())>1):     label_rgyrSD.append(t.RGyr.std())
        else:  label_rgyrSD.append(0)
        label_repl.append(np.unique(t.Repl.tolist(), return_counts=True))
        label_size.append(t.shape[0])

    
    starting_labels = df_ext[df_ext.time == 0].label.tolist() # list of clusters of all first frames in all trajectories
    starting_list, starting_count = np.unique(starting_labels, return_counts=True) # list of clusters that appear as first frame

    label2scale = label_rmsd # what value will be on x

    index_x_t = np.argsort(label2scale)
    # x and y positions of each cluster in the plot
    label_y = labels_list
    label_x = np.zeros((len(labels_list)),dtype = float)  

    max_x = max(label2scale)
    step_x = 1.0*max_x/max(label2scale)
    
    # order label_x: 
    # firt clasters that contain first frames of trajectories
    for i,s in enumerate(starting_list):
        label_x[s] = label2scale[s]*step_x
        label_y[s] = labels_list[s]
    # then the rest 
    j = 0
    for l in index_x_t:
        if (labels_list[l] not in starting_list):
            while (label_x[j] != 0):
                j += 1
                if(j == len(labels_list)): break
            if(j == len(labels_list)):break
            label_x[labels_list[l]] = label2scale[l]*step_x 
            label_y[labels_list[l]] = labels_list[l] #5+20*np.random.rand(1)[0]
        
    # s in label_x: print("x-position",s)    
    #color = ['r','b','forestgreen','lime','m','c','teal','orange','yellow','goldenrod','olive','tomato','salmon','seagreen']
    label_x = np.log10(label_x)
    x_tick_lable = []
    x_tick_pos = []
    for k in range(0,2):
        for ii,i in enumerate(range(pow(10,k),pow(10,k+1),pow(10,k))):  
            if(ii == 0): x_tick_lable.append(str(i))
            else: x_tick_lable.append("")
            x_tick_pos.append(np.log10(i))
            if(i > 25): break

               
    if "ligand" in df_ext.columns.tolist():
        if replica_show:
            lig_repl_list = np.unique(np.asarray(tuple(zip(df_ext.ligand,df_ext.Repl))),axis=0)
        else:
            lig_repl_list = np.unique(np.asarray(df_ext.ligand))
    else:
        lig_repl_list = np.unique(df_ext.Repl.tolist())
        
    #------------------------------------------    
    for i,r in  enumerate(lig_repl_list):# loop over ligand/replicas or ligand
        if "ligand" in df_ext.columns.tolist(): 
            if replica_show:
                rr0 = df_ext[df_ext.ligand == r[0]]  
                rr = rr0[rr0.Repl == r[1]]  
                name = r[0]+" "+r[1]
            else:
                rr = df_ext[df_ext.ligand == r]
                name = r
        else:   
            rr = df_ext[df_ext.Repl == r]
            name = r
        
        fig = plt.figure(figsize=(16, 5))
        gs = gridspec.GridSpec(1, 3, width_ratios=[6, 1, 1]) 
        ax = plt.subplot(gs[0])
        subset_labels = np.unique(rr.label.tolist())
        #  firs plot all cluster nodes
        ax.scatter(label_x[subset_labels],label_y[subset_labels],\
                   facecolors='none',color ='orange',marker='o',alpha=0.8,\
               s=2*np.asarray(label_size)[subset_labels])
        
        for txt in labels_list[subset_labels]:
            ax.annotate(txt, (label_x[txt],label_y[txt]),fontsize=12)
        
        for j,traj in  enumerate(np.unique(rr.Traj.tolist())): # loop over trajectories
            # select only trajectories where ligand completely dissociated
            rr_traj_r = rr[rr.Traj == traj]
            
            for repl  in np.unique(rr_traj_r.Repl.tolist()):
                rr_traj = rr_traj_r[rr_traj_r.Repl == repl] 
                if((np.asarray(rr_traj.RMSDl.tolist())[-1] > 10 ) \
                        and (np.asarray(rr_traj.RMSDl.tolist())[0] < 2.5)):
                    traj_list = np.asarray(rr_traj.label.tolist())  # list of the lables for all frame sequentially
                    traj_x = []
                    traj_y = []
                    for t in traj_list: # list of the cluster visited in the trajectory
                        traj_x.append(label_x[t])
                        traj_y.append(label_y[t])
                    color = "k"
                # show one selected trajectory (replica, trajectory)
                    lw = 0.5 # line width in the plot
                    color = 'k'
                    if cluster_show >= 0:
                        if (rr_traj.traj_cluster.tolist()[0] == cluster_show): 
                            lw ,color = 5,"dodgerblue"
                    ax.plot(traj_x,traj_y, color = color,linewidth=lw,alpha=0.5)
                # # let us plot final point in the trajectory
                    ax.scatter(traj_x[-1],traj_y[-1], color ='red',alpha=0.3,s=500)
                    ax.scatter(traj_x[0],traj_y[0], color ='green',alpha=0.3,s=500)
                
        ax.set_ylabel('Cluster', fontsize=20)
        ax.set_xlabel('log10(<RMSD>) /a.u.', fontsize=16) 
        plt.xticks(x_tick_pos,x_tick_lable, fontsize=16)
        ax.tick_params(labelsize=16)
        plt.title(name, fontsize=20)
        ax1 = plt.subplot(gs[1])
        ax1.errorbar(x=np.asarray(label_wat)[subset_labels],y=np.asarray(labels_list)[subset_labels],\
                 xerr=np.asarray(label_watSD)[subset_labels],color = 'k', alpha=0.6,marker='o',lw=0.5)
        ax1.set_xlabel('water #', fontsize=16) 
        ax1.tick_params(labelsize=16)
        ax2 = plt.subplot(gs[2])
        ax2.errorbar(x=np.asarray(label_rgyr)[subset_labels],y=np.asarray(labels_list)[subset_labels],\
                 xerr=np.asarray(label_rgyrSD)[subset_labels],color = 'k', alpha=0.6,marker='o',lw=0.5)
        ax2.set_xlabel('RGyr #', fontsize=16) 
        ax2.tick_params(labelsize=16)
#    plt.savefig(tr.PRJ_DIR+name+"-ifp.png")
    plt.show()
    return 
    """

   
    label_time = []
    label_timeSD = []
    label_rmsd = []
    label_rmsdSD = []
    label_wat = []
    label_watSD = []
    label_repl = []
    label_size = []
    label_rgyr = []
    label_rgyrSD = []
    label_length = []
    label_lengthSD = []


    labels_list = np.unique(df_ext.label.tolist())
    for l in labels_list:
        t = df_ext[df_ext.label == l]
        label_time.append(int(t.time.mean()))
        if(len(t.time.tolist())>1):     label_timeSD.append(int(t.time.std()))
        else:  label_timeSD.append(0)
        label_length.append(int(t.length.mean()))
        if(len(t.time.tolist())>1):     label_lengthSD.append(int(t.length.std()))
        else:  label_lengthSD.append(0)
        label_rmsd.append(t.RMSDl.mean())
        if(len(t.time.tolist())>1):     label_rmsdSD.append(t.RMSDl.std())
        else:  label_rmsdSD.append(0)
        label_wat.append(int(t.WAT.mean()))
        if(len(t.time.tolist())>1):     label_watSD.append(int(t.WAT.std()))
        else:  label_watSD.append(0)
        label_rgyr.append(t.RGyr.mean())
        if(len(t.time.tolist())>1):     label_rgyrSD.append(t.RGyr.std())
        else:  label_rgyrSD.append(0)
        label_repl.append(np.unique(t.Repl.tolist(), return_counts=True))
        label_size.append(t.shape[0])

    
    starting_labels = df_ext[df_ext.time == 0].label.tolist() # list of clusters of all first frames in all trajectories
    starting_list, starting_count = np.unique(starting_labels, return_counts=True) # list of clusters that appear as first frame

    label2scale = label_rmsd # what value will be on x

    index_x_t = np.argsort(label2scale)
    # x and y positions of each cluster in the plot
    label_y = labels_list
    label_x = np.zeros((len(labels_list)),dtype = float)  

    max_x = max(label2scale)
    step_x = 1.0*max_x/max(label2scale)
    
    # order label_x: 
    # firt clasters that contain first frames of trajectories
    for i,s in enumerate(starting_list):
        label_x[s] = label2scale[s]*step_x
        label_y[s] = labels_list[s]
    # then the rest 
    j = 0
    for l in index_x_t:
        if (labels_list[l] not in starting_list):
            while (label_x[j] != 0):
                j += 1
                if(j == len(labels_list)): break
            if(j == len(labels_list)):break
            label_x[labels_list[l]] = label2scale[l]*step_x 
            label_y[labels_list[l]] = labels_list[l] #5+20*np.random.rand(1)[0]
        
    # s in label_x: print("x-position",s)    
    #color = ['r','b','forestgreen','lime','m','c','teal','orange','yellow','goldenrod','olive','tomato','salmon','seagreen']
    label_x = np.log10(label_x)
    x_tick_lable = []
    x_tick_pos = []
    for k in range(0,2):
        for ii,i in enumerate(range(pow(10,k),pow(10,k+1),pow(10,k))):  
            if(ii == 0): x_tick_lable.append(str(i))
            else: x_tick_lable.append("")
            x_tick_pos.append(np.log10(i))
            if(i > 25): break

               
    if "ligand" in df_ext.columns.tolist():
        if replica_show:
            lig_repl_list = np.unique(np.asarray(tuple(zip(df_ext.ligand,df_ext.Repl))),axis=0)
        else:
            lig_repl_list = np.unique(np.asarray(df_ext.ligand))
    else:
        lig_repl_list = np.unique(df_ext.Repl.tolist())
        
    #------------------------------------------    
    for i,r in  enumerate(lig_repl_list):# loop over ligand/replicas or ligand
        if "ligand" in df_ext.columns.tolist(): 
            if replica_show:
                rr0 = df_ext[df_ext.ligand == r[0]]  
                rr = rr0[rr0.Repl == r[1]]  
                name = r[0]+" "+r[1]
            else:
                rr = df_ext[df_ext.ligand == r]
                name = r
        else:   
            rr = df_ext[df_ext.Repl == r]
            name = r
        
        fig = plt.figure(figsize=(16, 5))
        gs = gridspec.GridSpec(1, 4, width_ratios=[1, 6, 1, 1]) 
        ax = plt.subplot(gs[1])
        subset_labels = np.unique(rr.label.tolist())
        #  firs plot all cluster nodes
        ax.scatter(label_x[subset_labels],label_y[subset_labels],\
                   facecolors='none',color ='orange',marker='o',alpha=0.8,\
               s=2*np.asarray(label_size)[subset_labels])
        
        for txt in labels_list[subset_labels]:
            ax.annotate(txt, (label_x[txt],label_y[txt]),fontsize=18)
        
        for j,traj in  enumerate(np.unique(rr.Traj.tolist())): # loop over trajectories
            # select only trajectories where ligand completely dissociated
            rr_traj_r = rr[rr.Traj == traj]
            for repl  in np.unique(rr_traj_r.Repl.tolist()):
                rr_traj = rr_traj_r[rr_traj_r.Repl == repl] 
#                print("----",traj,np.asarray(rr_traj.RMSDl.tolist())[-1],\
#                                             np.asarray(rr_traj.RMSDl.tolist())[0])
                if((np.asarray(rr_traj.RMSDl.tolist())[-1] > 10 ) \
                        and (np.asarray(rr_traj.RMSDl.tolist())[0] < 5)):
                    traj_list = np.asarray(rr_traj.label.tolist())  # list of the lables for all frame sequentially
                    traj_x = []
                    traj_y = []
                    for t in traj_list: # list of the cluster visited in the trajectory
                        traj_x.append(label_x[t])
                        traj_y.append(label_y[t])
                    color = "k"
                # show one selected trajectory (replica, trajectory)
                    lw = 2 # line width in the plot
                    color = 'k'
                    if cluster_show >= 0:
                        if (rr_traj.traj_cluster.tolist()[0] == cluster_show): 
                            lw ,color = 5,"dodgerblue"
                    ax.plot(traj_x,traj_y, color = color,linewidth=lw,alpha=0.3)
                # # let us plot final point in the trajectory
                    ax.scatter(traj_x[-1],traj_y[-1], color ='red',alpha=0.3,s=500)
                    ax.scatter(traj_x[0],traj_y[0], color ='green',alpha=0.3,s=500)

                
        ax.set_ylabel('Cluster', fontsize=20)
        ax.set_xlabel('log10(<RMSD>) /a.u.', fontsize=20) 
        plt.xticks(x_tick_pos,x_tick_lable, fontsize=20)
        ax.tick_params(labelsize=20)
        plt.title(name, fontsize=20)
        plt.grid()
        
        ax1 = plt.subplot(gs[2])
        ax1.errorbar(x=np.asarray(label_wat)[subset_labels],y=np.asarray(labels_list)[subset_labels],\
                 xerr=np.asarray(label_watSD)[subset_labels],color = 'k', alpha=0.6,marker='o',lw=0.5)
        ax1.set_xlabel('water #', fontsize=20) 
        ax1.tick_params(labelsize=20)
        plt.grid()
        
        ax2 = plt.subplot(gs[3])
        ax2.errorbar(x=np.asarray(label_rgyr)[subset_labels],y=np.asarray(labels_list)[subset_labels],\
                 xerr=np.asarray(label_rgyrSD)[subset_labels],color = 'k', alpha=0.6,marker='o',lw=0.5)
        ax2.set_xlabel('RGyr #', fontsize=20) 
        ax2.tick_params(labelsize=20)
        plt.grid()
        
        ax0 = plt.subplot(gs[0])
        population = np.zeros((len(labels_list)))
        for lab in np.sort(np.unique(subset_labels)):
            population[lab] = rr[rr.label==lab].shape[0]
        ax0.errorbar(x=population/np.max(population),y=np.asarray(labels_list)\
                 ,color = 'k', alpha=0.6,marker='o',lw=1)
        ax0.set_xlabel('population', fontsize=20) 
        ax0.tick_params(labelsize=20)
        plt.grid()

#    plt.savefig(tr.PRJ_DIR+name+"-ifp.png")
    plt.show()
    return  (index_x_t)

######################################
#
#
#####################################
def Plot_cluster_IFP(df_ext,ligand,label_list,resi_name_list_sorted,ind_part_resi,delta_list=[0,0],chain_length = -1,file_save=""):  #[18,164]):



    df_ext_lig = df_ext[df_ext.ligand == ligand]

    contact_list = []
    population_list = []
    for label in label_list:
        df_ext_lig_lab = df_ext_lig[df_ext_lig.label == label]
        c = []
        p = []
        for l in df_ext_lig_lab.columns.tolist():
            try:
                if(chain_length > 0 and int(l[6:]) > chain_length): delta=delta_list[1]
                else: delta = delta_list[0]
                delta = 0
                new_list_col= l[3:6]+str(int(l[6:])+delta)
                if new_list_col in resi_name_list_sorted[ind_part_resi]:
                    if(df_ext_lig_lab[l].sum(axis=0)> 5):
 #                   print(l,new_list_col,df_ext_lig_lab[l].sum(axis=0))
                        c.append(new_list_col)
                        p.append(df_ext_lig_lab[l].sum(axis=0))
            except:
                pass
        contact_list.append(c)
        population_list.append(p)
    
    resi_list = []
    for l in contact_list:
        for c in l: resi_list.append(c)
        
    resi_list = np.unique(np.asarray(resi_list))

    n = []
    for t in resi_list:  n.append(t[3:])
    ind_sorted = np.argsort(np.asarray(n))  


    ar = np.zeros((len(label_list), len(resi_list)))
    for i,l in enumerate(label_list):
        for c,p in zip(contact_list[i],population_list[i]):
            ind = np.argwhere(resi_list[ind_sorted] == c)[0][0]
            ar[i,ind]=p
        
    fig = plt.figure(figsize = (16, 8),facecolor='w')
    ax = plt.subplot(1,1,1)
    plt.imshow(ar,cmap='Blues')
    plt.xticks(range(0,len(resi_list)),resi_list[ind_sorted],rotation=90,fontsize=16)
    plt.yticks(range(0,len(label_list)),label_list,rotation=90,fontsize=16)
    if file_save !="":plt.savefig(file_save,dpi=300)  
    else: plt.show()
#ax.set_yticklabels(label_list)
    return


########################################
#
#     Clusterign of trajectories based on the list of graph set for all trajectories
#
########################################
   
def cluster_trajectories(df_ext,labels_list, n_CL=10, RMSD_threshold=3,transition_matrix = False):
    """
    Parameters:
    df_ext - IFP database
    labels_list- practically a list of cluster numbers obtained from clastaring of all snapshots
    n_CL - number of clusters
    RMSD_threshold - RMSD threshold setting  which snapshot should be considered
    transition_matrix - defines either transition matrix between states (True) 
                        or just the frequency of visiting a state should be used  for clustering
    
    Returns:
    traj_all_clusters
    kmeansT - list of labels
    similarity_summary -  mean and std of pair-wise similarity in the cluster 
    """

    similarity_summary = []
    #  database may contain "ligand" column
    if "ligand" in df_ext.columns.tolist():
        lig_repl_list = np.unique(np.asarray(tuple(zip(df_ext.ligand,df_ext.Repl))),axis=0)
    else:
        lig_repl_list = np.unique(df_ext.Repl.tolist())
            
    graph = []  # will accomulate transition matrixes for all trajectories
    complete_name = []
    l = 0
    save_x = []  # in this array we will save exactly ligand, replica, and traj to be able to asign trajectory to a particular cluster
    
    for i,r in  enumerate(lig_repl_list):# loop over ligands and replicas
        if "ligand" in df_ext.columns.tolist(): 
            rr0 = df_ext[df_ext.ligand == r[0]]  # ligand sub-set
            rr = rr0[rr0.Repl == r[1]]           # ligand-replica sub-set
            name = r[0]+" "+r[1]
        else:   
            print("Error: missing column ligand in the database! ")  # just replica sub-set in the case if "ligand" column is absent

        for j,traj in  enumerate(np.unique(rr.Traj.tolist())): # loop over trajectories
            tt_all = rr[rr.Traj== traj] # trajectory sub-set
            tt = tt_all[tt_all.RMSDl > RMSD_threshold]  # part of the trajectory with RMSD above a certain threshold
            cluster, dens = np.unique(tt.label.tolist(), return_counts=True)  # clusters and their population
            if transition_matrix:
                #----------- transition matrix-----------------
                graph_traj = np.zeros((len(labels_list),len(labels_list)),dtype = float)  # transition matrix for one trajeectory
                for c,d in zip(cluster, dens):
                    graph_traj[c,c] = 1 #1.0*d/length   # cluster visited (diogonal terms of the transition matrix)
                for i,t in enumerate(tt.label.tolist()):
                    if(i > 0):
                        graph_traj[tt.label.tolist()[i],tt.label.tolist()[i-1]] = 1# += 1.0/length
                #------------------------------------------------
            else:
                #--------------- visiting matrix ------------------
                graph_traj = np.zeros((len(labels_list)),dtype = float) 
                for c,d in zip(cluster, dens):
                    graph_traj[c] = d/tt.shape[0]
                
            graph.append(graph_traj.flatten()) 
            complete_name.append(name+"_"+str(j))
            l += 1
            save_x.append(r[0]+" "+r[1]+" "+traj)
    
    X = np.asarray(graph)
    
    print("===== Matrix for clustering ====")
    print("that containes flattered graphs matrixes for all RAMD trajectories (shape: X.shape)")
    print("Database:",df_ext.shape,"Matrix",X.shape,"altougether trajectoies:",l)
    
    similarity_matrix = np.zeros((X.shape[0],X.shape[0]),dtype = float)
    # compute and plot similarity matrix
    for i in range(0,X.shape[0]):
        sim = []
        for j in range(i+1,X.shape[0]):
            similarity_matrix[i,j] = np.linalg.norm(graph[i]-graph[j])#np.sum(np.multiply(graph[i],graph[j]))
            similarity_matrix[j,i] = similarity_matrix[i,j] 

            

    kmeansT = KMeans(n_clusters=n_CL, random_state=0, n_init = 20, max_iter = 500).fit(X)
    
    Y = []
    for l in np.unique(kmeansT.labels_):
        x0 = []
        for i in range(0,X.shape[0]):
            if(kmeansT.labels_[i] == l): x0.append(X[i])  
        Y.append(np.mean(np.asarray(x0),axis = 0))
    plt.figure(figsize = (10,10))
    g = sns.heatmap(np.asarray(Y).T,linewidth = 0, cmap ="YlGnBu")
    plt.show()
    
    # now we will make a list that asigns labels to each snapshot in the database
    traj_clust_labels = np.zeros((df_ext.shape[0]),dtype = np.int8)
    save_x = np.asarray(save_x)
    j = 0
    for i in range(0,df_ext.shape[0]):
        l_r_t = df_ext.iloc[i]["ligand"]+" "+df_ext.iloc[i]["Repl"]+" "+df_ext.iloc[i]["Traj"]
        if(l_r_t in save_x):
            traj_clust_labels[i] = kmeansT.labels_[np.argwhere(save_x == l_r_t)]
    
    # Visualize clusters

    mask = kmeansT.labels_
    sim_p_w_m_mean = []
    sim_p_w_m_std = []  
    sim_p_w_d_mean = []
    sim_p_w_d_std = []  
    for l in np.unique(kmeansT.labels_):
            index = np.argwhere(kmeansT.labels_ == l).flatten()
            indexF = np.argwhere(kmeansT.labels_ != l).flatten()
            m = similarity_matrix[index].T[index].T
            d = similarity_matrix[index].T[indexF].T
            sim_p_w_m_mean.append(np.round(np.mean(m.flatten()),1))
            sim_p_w_m_std.append(np.round(np.std(m.flatten()),1))
            sim_p_w_d_mean.append(np.round(np.mean(d.flatten()),1))
            sim_p_w_d_std.append(np.round(np.std(d.flatten()),1))
            similarity_summary.append( (np.round(np.mean(m.flatten()),1),np.round(np.std(m.flatten()),1)))
            md = np.concatenate((m,d),axis=1)
    fig = plt.figure(figsize=(16, 8))
    gs = gridspec.GridSpec(2, 2,hspace=0.3,height_ratios=[1,1],width_ratios=[8,1]) 
    ax1 = plt.subplot(gs[0])
    ax1.errorbar(x=np.unique(kmeansT.labels_),y=sim_p_w_m_mean,yerr=sim_p_w_m_std,color ="red", \
             alpha=0.6,marker='o',lw=0.5,ms = 3, label = "in cluster",linestyle = "")
    ax1.errorbar(x=np.unique(kmeansT.labels_),y=sim_p_w_d_mean,yerr=sim_p_w_d_std,color ="blue", \
             alpha=0.6,marker='o',lw=0.5,ms = 3, label = "off",linestyle = "")
    ax1.tick_params(labelsize=16)
    ax1.legend()
    plt.show()
        

    l = 0
    all_traj_labels = []
    
    for i,r in  enumerate(lig_repl_list):# loop over replica
        if "ligand" in df_ext.columns.tolist(): 
            rr0 = df_ext[df_ext.ligand == r[0]]  
            rr = rr0[rr0.Repl == r[1]]  
            name = r[0]+" "+r[1]
        else:   
            rr = df_ext[df_ext.Repl == r]
            name = r
        traj_labels = []
        for j,traj in  enumerate(np.unique(rr.Traj.tolist())): # loop over trajectories
            for frame in range(0,len(rr[rr.Traj == traj].Traj.tolist())): 
                all_traj_labels.append(kmeansT.labels_[l])
                traj_labels.append(kmeansT.labels_[l])
            l += 1
 #       print(name,np.unique(np.asarray(traj_labels))
    
    df_ext["traj_cluster"] = np.asarray(all_traj_labels)

    return (traj_clust_labels,similarity_summary)   

########################################################################
# Plotting trajectories 
########################################################################

def plot_graph_test(df_ext,cluster_show=-1,replica_show = False,file_save = ""):
    """
    Parameters:
    df_ext - IFP database
    cluster_show - cluster of trajectories to be highlighed 
    Returns:
    """
   
    label_time = []
    label_timeSD = []
    label_rmsd = []
    label_rmsdSD = []
    label_wat = []
    label_watSD = []
    label_repl = []
    label_size = []
    label_rgyr = []
    label_rgyrSD = []
    label_length = []
    label_lengthSD = []


    labels_list = np.unique(df_ext.label.tolist())
    for l in labels_list:
        t = df_ext[df_ext.label == l]
        label_time.append(int(t.time.mean()))
        if(len(t.time.tolist())>1):     label_timeSD.append(int(t.time.std()))
        else:  label_timeSD.append(0)
        label_length.append(int(t.length.mean()))
        if(len(t.time.tolist())>1):     label_lengthSD.append(int(t.length.std()))
        else:  label_lengthSD.append(0)
        label_rmsd.append(t.RMSDl.mean())
        if(len(t.time.tolist())>1):     label_rmsdSD.append(t.RMSDl.std())
        else:  label_rmsdSD.append(0)
        label_wat.append(int(t.WAT.mean()))
        if(len(t.time.tolist())>1):     label_watSD.append(int(t.WAT.std()))
        else:  label_watSD.append(0)
        label_rgyr.append(t.RGyr.mean())
        if(len(t.time.tolist())>1):     label_rgyrSD.append(t.RGyr.std())
        else:  label_rgyrSD.append(0)
        label_repl.append(np.unique(t.Repl.tolist(), return_counts=True))
        label_size.append(t.shape[0])

    
    starting_labels = df_ext[df_ext.time == 0].label.tolist() # list of clusters of all first frames in all trajectories
    starting_list, starting_count = np.unique(starting_labels, return_counts=True) # list of clusters that appear as first frame

    label2scale = label_rmsd # what value will be on x

    index_x_t = np.argsort(label2scale)
    # x and y positions of each cluster in the plot
    label_y = labels_list
    label_x = np.zeros((len(labels_list)),dtype = float)  

    max_x = max(label2scale)
    step_x = 1.0*max_x/max(label2scale)
    
    # order label_x: 
    # firt clasters that contain first frames of trajectories
    for i,s in enumerate(starting_list):
        label_x[s] = label2scale[s]*step_x
        label_y[s] = labels_list[s]
    # then the rest 
    j = 0
    for l in index_x_t:
        if (labels_list[l] not in starting_list):
            while (label_x[j] != 0):
                j += 1
                if(j == len(labels_list)): break
            if(j == len(labels_list)):break
            label_x[labels_list[l]] = label2scale[l]*step_x 
            label_y[labels_list[l]] = labels_list[l] #5+20*np.random.rand(1)[0]
        
    # s in label_x: print("x-position",s)    
    #color = ['r','b','forestgreen','lime','m','c','teal','orange','yellow','goldenrod','olive','tomato','salmon','seagreen']
    label_x = np.log10(label_x)
    x_tick_lable = []
    x_tick_pos = []
    for k in range(0,2):
        for ii,i in enumerate(range(pow(10,k),pow(10,k+1),pow(10,k))):  
            if(ii == 0): x_tick_lable.append(str(i))
            else: x_tick_lable.append("")
            x_tick_pos.append(np.log10(i))
            if(i > 25): break

               
    if "ligand" in df_ext.columns.tolist():
        if replica_show:
            lig_repl_list = np.unique(np.asarray(tuple(zip(df_ext.ligand,df_ext.Repl))),axis=0)
        else:
            lig_repl_list = np.unique(np.asarray(df_ext.ligand))
    else:
        lig_repl_list = np.unique(df_ext.Repl.tolist())
        
    #------------------------------------------    
    for i,r in  enumerate(lig_repl_list):# loop over ligand/replicas or ligand
        if "ligand" in df_ext.columns.tolist(): 
            if replica_show:
                rr0 = df_ext[df_ext.ligand == r[0]]  
                rr = rr0[rr0.Repl == r[1]]  
                name = r[0]+" "+r[1]
            else:
                rr = df_ext[df_ext.ligand == r]
                name = r
        else:   
            rr = df_ext[df_ext.Repl == r]
            name = r
        
        fig = plt.figure(figsize=(16, 5))
        gs = gridspec.GridSpec(1, 4, width_ratios=[1, 6, 1, 1]) 
        ax = plt.subplot(gs[1])
        subset_labels = np.unique(rr.label.tolist())
        #  firs plot all cluster nodes
        ax.scatter(label_x[subset_labels],label_y[subset_labels],\
                   facecolors='none',color ='orange',marker='o',alpha=0.8,\
               s=2*np.asarray(label_size)[subset_labels])
        
        for txt in labels_list[subset_labels]:
            ax.annotate(txt, (label_x[txt],label_y[txt]),fontsize=18)
        
        for j,traj in  enumerate(np.unique(rr.Traj.tolist())): # loop over trajectories
            # select only trajectories where ligand completely dissociated
            rr_traj_r = rr[rr.Traj == traj]
            for repl  in np.unique(rr_traj_r.Repl.tolist()):
                rr_traj = rr_traj_r[rr_traj_r.Repl == repl] 
#                print("----",traj,np.asarray(rr_traj.RMSDl.tolist())[-1],\
#                                             np.asarray(rr_traj.RMSDl.tolist())[0])
                if((np.asarray(rr_traj.RMSDl.tolist())[-1] > 10 ) \
                        and (np.asarray(rr_traj.RMSDl.tolist())[0] < 5)):
                    traj_list = np.asarray(rr_traj.label.tolist())  # list of the lables for all frame sequentially
                    traj_x = []
                    traj_y = []
                    for t in traj_list: # list of the cluster visited in the trajectory
                        traj_x.append(label_x[t])
                        traj_y.append(label_y[t])
                    color = "k"
                # show one selected trajectory (replica, trajectory)
                    lw = 0.5 # line width in the plot
                    color = 'k'
                    if cluster_show >= 0:
                        if (rr_traj.traj_cluster.tolist()[0] == cluster_show): 
                            lw ,color = 5,"dodgerblue"
                    ax.plot(traj_x,traj_y, color = color,linewidth=lw,alpha=0.5)
                # # let us plot final point in the trajectory
                    ax.scatter(traj_x[-1],traj_y[-1], color ='red',alpha=0.3,s=500)
                    ax.scatter(traj_x[0],traj_y[0], color ='green',alpha=0.3,s=500)

                
        ax.set_ylabel('Cluster', fontsize=18)
        ax.set_xlabel('log10(<RMSD>) /a.u.', fontsize=18) 
        plt.xticks(x_tick_pos,x_tick_lable, fontsize=18)
        ax.tick_params(labelsize=18)
        plt.title(name, fontsize=18)
        plt.grid()
        
        ax1 = plt.subplot(gs[2])
        ax1.errorbar(x=np.asarray(label_wat)[subset_labels],y=np.asarray(labels_list)[subset_labels],\
                 xerr=np.asarray(label_watSD)[subset_labels],color = 'k', alpha=0.6,marker='o',lw=0.5)
        ax1.set_xlabel('water #', fontsize=18) 
        ax1.tick_params(labelsize=18)
        plt.grid()
        
        ax2 = plt.subplot(gs[3])
        ax2.errorbar(x=np.asarray(label_rgyr)[subset_labels],y=np.asarray(labels_list)[subset_labels],\
                 xerr=np.asarray(label_rgyrSD)[subset_labels],color = 'k', alpha=0.6,marker='o',lw=0.5)
        ax2.set_xlabel('RGyr #', fontsize=18) 
        ax2.tick_params(labelsize=18)
        plt.grid()
        
        ax0 = plt.subplot(gs[0])
        population = np.zeros((len(labels_list)))
        for lab in np.sort(np.unique(subset_labels)):
            population[lab] = rr[rr.label==lab].shape[0]
        ax0.errorbar(x=population/np.max(population),y=np.asarray(labels_list)\
                 ,color = 'k', alpha=0.6,marker='o',lw=1)
        ax0.set_xlabel('population', fontsize=18) 
        ax0.tick_params(labelsize=18)
        plt.grid()

#    plt.savefig(tr.PRJ_DIR+name+"-ifp.png")
    if file_save != "": plt.savefig(file_save,dpi=300)  
    else:    plt.show()
    return  (index_x_t)




########################################################################
# Sorting of residue by number
########################################################################

def get_resn_list(li,lab):
    ret = []
    n = []
    for i in li:
        if i[:2] == lab:
            ret.append(i)
            n.append(int(i[6:]))
    ind_sorted = np.argsort(np.asarray(n))
    return np.asarray(ret)[ind_sorted]


########################################################################
# Program for reading IFP databases
# Additionally column with ligand name is added
# and COM column is splitted to COM_x, COM_y, COM_z
########################################################################
def standard_IFP(unpickled_dfi,ligandsi):
    """
    Parameters:
    dictionary of files with ITP databases {name1:file_path1[,name2:filepath2],...}
    Returns:
    combined IFP database
    """
#    unpickled_dfi = []
#    ligandsi = []
#    for lig in list_IFP:
#        unpickled_dfi.append(pd.read_pickle(list_IFP[lig]))
#        ligandsi.append(lig)

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

########################################################################
########################################################################
def separate_IFP(complete_list_IFP):
    resi_list = []
    ifp_list = []
    resi_name_list = []
    for i,name in enumerate(complete_list_IFP):
        if(name[2]== "_"):
            resi = name[6:]
            if resi not in resi_list: 
                resi_list.append(resi)
                resi_name_list.append(name[3:])
                ifp_list.append([0,0,0,0,0])
            ind = np.argwhere(np.asarray(resi_list) == resi)[0][0]
            if(name[0:2] == "AR"):  ifp_list[ind][0]=1
            if(name[0:2] == "HY"):  ifp_list[ind][1]=1
            if(name[0:2] == "HD" or name[0:2] == "HA"):  ifp_list[ind][2] = 1
            if(name[0:2] == "WB"):  ifp_list[ind][3] = 1
            if(name[0:2] == "RE"):  ifp_list[ind][4] = 1
        else:
            print(name, "- skip no-IFP property")
    ind_sorted = np.argsort(np.asarray(resi_list).astype(int))
    resi_list_sorted = np.sort(np.asarray(resi_list).astype(int)).astype(str)
    resi_list_sorted = np.asarray(resi_list)[ind_sorted]  
    resi_name_list_sorted = np.asarray(resi_name_list)[ind_sorted]
    return(resi_list_sorted,resi_name_list_sorted,ifp_list)

########################################################################
########################################################################
def get_from_prop(list_x, df,list_l= [],threshold = 0.1):
    """
    """
    if len(list_l) == 0:
        list_l = np.unique(df.ligand.tolist())
    ar = []
    ar_SD = []
    for ligand in np.unique(df.ligand.tolist()):
        df_ligand = df[df.ligand == ligand]
        if ligand in list_l:
            ar_repl = []
            for Repl in np.unique(df_ligand.Repl.tolist()):
                df_ligand_Repl = df_ligand[df_ligand.Repl == Repl]
                repl_mean = df_ligand_Repl[list_x].mean().values
                ar_repl.append(repl_mean)
            ar.append(np.mean(np.asarray(ar_repl),axis=0))
            ar_SD.append(np.std(np.asarray(ar_repl),axis=0))
    ar= np.asarray(ar)
    ar_SD= np.asarray(ar_SD)
    ind = np.where(np.mean(ar,axis=0)<threshold)
    ar=np.delete(ar,ind,1)
    ar_SD=np.delete(ar_SD,ind,1)
    x = np.delete(list_x,ind)
    return(ar,ar_SD,x)

########################################################################
########################################################################
def unify_resi(list_resi, df,resi_list_sorted,list_l= [], threshold=3):
    """
    list_resi - a complete list of IFP contacts to be considered
    resi_list_sorted - sorted residue numbers to be included in the IFP matrix
    list_l - list of ligands to be considered
    """
    if len(list_l) == 0:
        list_l = np.unique(df.ligand.tolist())
    ar = []
    ar_SD = []
    for ligand in list_l:
        df_ligand = df[df.ligand == str(ligand)]
        comx = df_ligand[df_ligand.time == 0].COM_x.mean(axis=0)
        comy = df_ligand[df_ligand.time == 0].COM_y.mean(axis=0)
        comz = df_ligand[df_ligand.time == 0].COM_z.mean(axis=0)
#        print(">>>>>>>>>>>>>>>>>>>>>>",comx,comy,comz)
        t = (df_ligand.COM_x-comx)*(df_ligand.COM_x-comx)+\
        (df_ligand.COM_y-comy)*(df_ligand.COM_y-comy)+(df_ligand.COM_z-comz)*(df_ligand.COM_z-comz)  
        if(threshold > 0):
            df_ligand_diss = df_ligand[t > threshold*threshold]
        else:
            df_ligand_diss = df_ligand[t < threshold*threshold]
        ar_repl = []
        """
        for Repl in np.unique(df_ligand.Repl.tolist()):
            df_ligand_Repl = df_ligand[df_ligand.Repl == Repl]
            repl_mean = df_ligand_Repl[list_resi].mean().values
            ar_repl.append(repl_mean)
        ar.append(np.mean(np.asarray(ar_repl),axis=0))
        ar_SD.append(np.std(np.asarray(ar_repl),axis=0))
        """
        ar.append(np.asarray(df_ligand_diss[list_resi].mean().values))
        ar_SD.append(np.asarray(df_ligand_diss[list_resi].std().values))
    ar= np.asarray(ar)
    ar_SD= np.asarray(ar_SD)
    x = list_resi
    """
    ind = np.where(np.mean(ar,axis=0)<threshold)
    if len(ind) > 9:
        ar=np.delete(ar,ind,1)
        ar_SD=np.delete(ar_SD,ind,1)
        x = np.delete(list_resi,ind)
    """
        
    ar_complete = np.zeros((len(list_l),len(resi_list_sorted)),dtype = float)
    ar_SD_complete = np.zeros((len(list_l),len(resi_list_sorted)),dtype = float)
    for k,ligand in enumerate(list_l):
        for i,xx in enumerate(x):
            try:
                ind = np.argwhere(xx[6:] == resi_list_sorted)
                ar_complete[k][ind] = ar[k][i]
                ar_SD_complete[k][ind] = ar_SD[k][i]
            except:
                pass  # this is in the case if we  left out part of residues in resi_list_sorted
    return(ar_complete,ar_SD_complete)

########################################################################
########################################################################
def ar_complete_ligand(ligand,df_tot,resi_list_sorted,properties=["RE","AR","HD","HA","HY","WB"]):
    """
    """
    df_ligand = df_tot[df_tot.ligand == ligand]
    ar_complete = np.zeros((len(properties),len(resi_list_sorted)),dtype = float)
    ar_SD_complete = np.zeros((len(properties),len(resi_list_sorted)),dtype = float)
    
    for k,pr in enumerate(properties):
        list_x = get_resn_list(df_ligand.columns.tolist(),pr)
        ar_repl = []
        for Repl in np.unique(df_ligand.Repl.tolist()):
            df_ligand_Repl = df_ligand[df_ligand.Repl == Repl]
            repl_mean = df_ligand_Repl[list_x].mean().values
            ar_repl.append(repl_mean)
        ar= np.mean(np.asarray(ar_repl),axis=0)
        ar_SD = (np.std(np.asarray(ar_repl),axis=0))
        for i,xx in enumerate(list_x):
            ind = np.argwhere(xx[6:] == resi_list_sorted)
            ar_complete[k][ind] = ar[i]
            ar_SD_complete[k][ind] = ar_SD[i]
            
    return(ar_complete,ar_SD_complete)


########################################################################
########################################################################
def read_databases(d,name_template,delta_list=[0,0],res_ren=-1):
    unpickled_dfi = []
    ligandsi = []
    list_IFP = {}    
    new_list_col = []
    for i,lig_pkl in enumerate(glob.glob(d+name_template)):
        name= lig_pkl[len(d):len(d)+8]
        print(lig_pkl,name)
        list_IFP.update( {name:lig_pkl})
        df_lig= pd.read_pickle(lig_pkl)
        list_col = df_lig.columns.tolist()
        new_list_col = [] # we need this array only once
        #--- to re-number residues
        for l in list_col:
            if(l[2] == "_"):  
                if((res_ren > 0) and (int(l[res_ren:]) > 200)): delta=delta_list[1]
                else: delta = delta_list[0]
                new_list_col.append(l[:6]+str(int(l[6:])+delta))
            else: new_list_col.append(l)
        df_lig.columns = new_list_col
        unpickled_dfi.append( df_lig)
        ligandsi.append(name)
    if len(new_list_col) <= 0:
        print("There is no files in :",d+name_template)
    df_tot = standard_IFP(unpickled_dfi,ligandsi)

    return(df_tot,ligandsi,new_list_col)


########################################################################
########################################################################
def clean_ramd(df_tot,threshold = 0.9):
    df_tot_new = pd.DataFrame(columns=df_tot.columns.tolist())
    for ligand in np.unique(df_tot.ligand):
        df_tot_ligand= df_tot[df_tot.ligand == ligand]
        for Repl in np.unique(df_tot_ligand.Repl):
            df_tot_ligand_Repl = df_tot_ligand[df_tot_ligand.Repl == Repl]
            list_out = 0
            for Traj in np.unique(df_tot_ligand_Repl.Traj):
                df_tot_ligand_Repl_Traj = df_tot_ligand_Repl[df_tot_ligand_Repl.Traj == Traj]
                if((df_tot_ligand_Repl_Traj.COM_z.tolist()[0])*threshold > df_tot_ligand_Repl_Traj.COM_z.tolist()[-1]):
                    print("Dropped",ligand, Repl, Traj)
                    list_out += 1
                else:
                    df_tot_new = pd.concat([df_tot_new,df_tot_ligand_Repl_Traj])
            if(list_out > 5): print(ligand,":   ",Repl)
    print(df_tot.shape, df_tot_new.shape)
    return(df_tot_new)



##########################################################################
######################################################################
def GRID_PRINT(file_name,pdrv,gr_orgn,gr_dim,grid_stp):
    header = "#  density  \n"
    header += "object 1 class gridpositions counts %3i" %gr_dim[0]+" %3i" %gr_dim[1]+" %3i" %gr_dim[2]+"\n"
    header += "origin %5.3f" %gr_orgn[0]+" %5.3f" %gr_orgn[1]+" %5.3f" %gr_orgn[2]+"\n"
    header += "delta %5.2f" %(grid_stp)+" 0 0 \n"
    header += "delta 0 %5.2f" %(grid_stp)+" 0 \n"
    header += "delta 0 0 %5.2f" %(grid_stp)+" \n"
    header += "object 2 class gridconnections counts %5i" %(gr_dim[0])+" %5i" %(gr_dim[1])+" %5i" %(gr_dim[2])+"\n"
    header += "object 3 class array type double rank 0 items %7i" %(gr_dim[0]*gr_dim[1]*gr_dim[2])+" data follows \n"

    check_range = int((3.0/grid_stp) + 1)

    output = []
    count = 0
    for i in pdrv.reshape(-1):
        output.append("%12.3e" %(i))
        count += 1
        if count%3 == 0:
            output.append("\n")
            
    with open(file_name,"w") as f:
        f.write(header)
        f.write("".join(output))
    return


##################################
################################
def Map_3D_grid(df_tot_to_save,filename):
    COM_x = []
    COM_y = []
    COM_z = []
    for x in df_tot_to_save.COM.values:
        COM_x.append(x[0])
        COM_y.append(x[1])
        COM_z.append(x[2])
    COM_x = np.asarray(COM_x)
    COM_y = np.asarray(COM_y)
    COM_z = np.asarray(COM_z)
    grid_mm_x = [COM_x.min(),COM_x.max()]
    grid_mm_y = [COM_y.min(),COM_y.max()]
    grid_mm_z = [COM_z.min(),COM_z.max()]
    grid_step = 1
    grid_dim= [int((grid_mm_x[1]-grid_mm_x[0])/grid_step+1),int((grid_mm_y[1]-grid_mm_y[0])/grid_step+1),int((grid_mm_z[1]-grid_mm_z[0])/grid_step+1)]
    grid = np.zeros((grid_dim),dtype=float)
    for (x,y,z) in zip(COM_x,COM_y,COM_z):
        ix= int((x-COM_x.min())/grid_step)
        iy= int((y-COM_y.min())/grid_step)
        iz= int((z-COM_z.min())/grid_step)
        grid[ix,iy,iz] += 1
    grid_origin = [grid_mm_x[0],grid_mm_y[0],grid_mm_z[0]]    
    GRID_PRINT(filename,grid,grid_origin,grid_dim,grid_step)
    return


##################################
################################


from matplotlib.patches import ArrowStyle
from matplotlib.patches import Ellipse
from scipy.spatial import distance
def plot_graph_New(df_ext,file_save = "",draw_round = False,merge=True):
    """
    Parameters:
    df_ext - IFP database
    Returns:
    """
   
    label_rmsd = []    # rmsd
    label_com = []    # COM
    label_size = []


    labels_list,nodes = np.unique(df_ext.label.values,return_counts= True)
    eidges = np.zeros((labels_list.shape[0],labels_list.shape[0]),dtype = float)
    coms = np.zeros((labels_list.shape[0],labels_list.shape[0]),dtype = float)

    
    for i,l in enumerate(labels_list):
        t = df_ext[df_ext.label == l]
        label_rmsd.append(t.RMSDl.mean())
        label_com.append(np.array((t.COM_x.mean(),t.COM_y.mean(),t.COM_z.mean())))
        label_size.append(t.shape[0])
        for j in range(0,i):
            coms[i,j] = distance.euclidean(label_com[i],label_com[j])
    
    for l,(df_label,df_time) in enumerate(zip(df_ext.label.values,df_ext.time.values)):
        if df_time != 0: 
            if(df_ext.label.values[l-1] != df_label):
                eidges[df_ext.label.values[l-1],df_label] += labels_list.shape[0]/df_ext.label.values.shape[0] 
    
            
 #   print(np.max(eidges), eidges[eidges > 0.5*np.max(eidges)])
    indx_first_com = np.argwhere((np.asarray(label_rmsd) == min(label_rmsd)))[0][0]
    dist_com = []
    for i,l in enumerate(labels_list): dist_com.append(np.round(distance.euclidean(label_com[i],label_com[indx_first_com]),2))
    dist_com = (10*np.asarray(dist_com)/np.max(dist_com)).astype(int)
    
    print(indx_first_com, min(label_rmsd),dist_com)
    fig = plt.figure(figsize = (6,2),facecolor='w') 
    gs = gridspec.GridSpec(1,2, width_ratios=[ 1,1],wspace=0.08) 
    ax = plt.subplot(gs[0]) 
    plt.title("Transition density")
    plt.imshow(eidges,cmap='Blues')
    ax = plt.subplot(gs[1])
    plt.title("Flow")
    flow = eidges-eidges.T
    plt.imshow(flow,cmap='Reds')
    plt.plot()


    starting_labels = df_ext[df_ext.time == 0].label.values # list of clusters of all first frames in all trajectories
    starting_list, starting_count = np.unique(starting_labels, return_counts=True) # list of clusters that appear as first frame

    #------------ oeder cluater position-------------

    label2scale = label_rmsd # what value will be on x
    label2order = labels_list #nodes #np.roll(labels_list,1) #nodes # what value will be on y

    index_x_t = np.argsort(label2scale)
    # x and y positions of each cluster in the plot
    label_x = np.zeros((len(labels_list)),dtype = float)  
    label_y = np.zeros((len(labels_list)),dtype = float)  
    color_com = np.zeros((len(labels_list)),dtype = float)  

    
    # order label_x: 
    max_x = max(label2scale)
    step_x = 1.0*max_x/max(label2scale)
    # firt clasters that contain first frames of trajectories
    for i,s in enumerate(starting_list):
        label_x[s] = label2scale[s]*step_x
        label_y[s] = label2order[s]
        color_com[s] = dist_com[s]
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
            color_com[labels_list[l]] = dist_com[l]
             
    # set logarythmic scale
    label_x = np.log10(label_x)
    x_tick_lable = []
    x_tick_pos = []
    for k in range(0,2):
        for ii,i in enumerate(range(pow(10,k),pow(10,k+1),pow(10,k))):  
 #           if(ii == 0): x_tick_lable.append(str(i))
 #           else: x_tick_lable.append("")
            x_tick_lable.append(str(i))
            x_tick_pos.append(np.log10(i))
            if(i > 25): break
      
    
    if draw_round:
        alpha = 0.9*2*3.14*label_x/np.max(label_x)
        alpha_regular = 0.9*2*3.14*np.asarray(x_tick_pos)/max(x_tick_pos)
        label_y = np.sin(alpha)
        label_x = np.cos(alpha)
        fig = plt.figure(figsize=(10, 10))
        gs = gridspec.GridSpec(1, 1) #, width_ratios=[1, 1]) 
        ax = plt.subplot(gs[0])
        plt.scatter(x=np.cos(alpha_regular),y=np.sin(alpha_regular), c='k',s=10)
        for l,p in zip(x_tick_lable,x_tick_pos):
            ax.annotate(str(l)+"A", (1.2*np.cos(0.9*2*3.14*p/max(x_tick_pos)),np.sin(0.9*2*3.14*p/max(x_tick_pos))),fontsize=14,color="gray")
        plt.xlim(-1.3,1.3)
        plt.ylim(-1.3,1.3)
    else:
        fig = plt.figure(figsize=(14, 10))
        gs = gridspec.GridSpec(1, 1) #, width_ratios=[1, 1]) 
        ax = plt.subplot(gs[0])
        ax.set_ylabel('Cluster', fontsize=18)
        ax.set_xlabel('<RMSD> /Angstrom', fontsize=18) 
        plt.xticks(x_tick_pos,x_tick_lable, fontsize=18)
        ax.tick_params(labelsize=18)
        plt.grid()

    
    el = Ellipse((2, -1), 0.4, 0.4)    
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
            if (eidges[l,n] > 0) :
                if  (np.abs((label_rmsd[l] - label_rmsd[n])) > 0.5* min(label_rmsd)) or (draw_round == False):
                    ax.annotate("", xy=xy, xycoords='data',
                        xytext=xytext, textcoords='data',
                        size=eidges[l,n]*500,
                        arrowprops=dict(arrowstyle="Fancy,head_length=0.2, head_width=0.4, tail_width=0.2", 
                                fc="orange", ec="none", alpha=0.2 ,
                                connectionstyle="arc3,rad=-0.5"),
                        )
                if  (np.abs((label_rmsd[l] - label_rmsd[n])) > 0.5* min(label_rmsd)) or (draw_round == False):
                    ax.annotate("", xy=xytext, xycoords='data',
                        xytext=xy, textcoords='data',
                        size=eidges[l,n]*500,
                        arrowprops=dict(arrowstyle="Fancy,head_length=0.2, head_width=0.4, tail_width=0.2", 
                                fc="orange", ec="none", alpha=0.2 ,
                                connectionstyle="arc3,rad=-0.5"),
                        )
            #  the flow
            flow = eidges[l,n] - eidges[n,l]  # flow l ----> n
            if (flow > 0) :
                a = l
                b = n
            else:
                a = n
                b = l                
            xy=(label_x[b],label_y[b])
            xytext=(label_x[a],label_y[a])   
            ax.annotate("", xy=xy, xycoords='data',
                        xytext=xytext, textcoords='data',
                        size=np.abs(flow)*4000,
                        arrowprops=dict(arrowstyle="Simple,head_length=0.2, head_width=0.4, tail_width=0.2", 
                                fc="0.6", ec="none", alpha=0.8 ,
                                connectionstyle="arc3,rad=-0.5"),
                        )

    for i,txt in enumerate(labels_list):
            ax.annotate(txt, (label_x[txt],label_y[txt]+0.05*pow(i,0.5)),fontsize=18)
            
    ax.scatter(label_x,label_y,facecolors='none',c =color_com,edgecolors="k", s=2*np.asarray(label_size), cmap='Oranges')
    
    if file_save != "": plt.savefig(file_save,dpi=300)  
    else:    plt.show()
    return
