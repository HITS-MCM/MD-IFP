import glob, os
import sys
import numpy as np

import pandas as pd

from matplotlib import *
from matplotlib import cm
import matplotlib.ticker
import  pylab as plt

def rank_IFP_resi(df,ifp_type=['AR','HY','HA','HD','HL','IP','IN',"IO","WB"]):
    """    
    this function extracts and ranks by the residue number IFP list  from the IFP table     
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



def Plot_IF_trajectory(df_tot,ifp_type = np.asarray(['AR','HA','HD','HY','IP','IN',"IO","WB"]),head_tail=-1,save_file = ""):
    """
    plotting IFPs: (i) average over complete trajectory, (ii) average over N fiest frames, (iii) averaged over last N frames
    Parameters:
    head_tail = N - number of firt/last frames to comute average
    ifp_type - list of IFP types to be considered
    df_tot - IFP pamndas frame
    """
    
    columns_IFP,columns_RE= rank_IFP_resi(df_tot, ifp_type)
    columns_sel = columns_IFP
    if head_tail < 0:
        head_tail = int(df_tot.shape[0]/3)       
    X = []
    df = df_tot[columns_sel] 
  #  columns_sel = columns_sel[df.mean().values > 0.01]
  #  columns_HB=[]
  #  [columns_HB.append(c) for c in columns_sel if (c[:2] == "HD" or c[:2] == "HA")]
  #  columns_HB = np.asarray(columns_HB)    

    fig = plt.figure(figsize=(14,3),dpi=150)

    df = df_tot[np.append(columns_sel,"time")] 
  #  n_hb = len(columns_HB[(df[columns_HB].mean()> 0.75).values])
    plt.bar(range(0,len(columns_sel)),df[df.time < head_tail][columns_sel].mean(),alpha=0.6,label="first "+str(head_tail)+" frames")
    plt.bar(range(0,len(columns_sel)),df[df.time > df.shape[0]-head_tail][columns_sel].mean(),alpha=0.6,label="last "+str(head_tail)+" frames")
    plt.bar(np.asarray(range(0,len(columns_sel))),df[columns_sel].mean(),color="",label="all frames",edgecolor ='k',hatch="/")
    plt.xticks(range(0,len(columns_sel)), columns_sel, rotation='vertical',fontsize=10)
    plt.legend(fontsize=10, loc = 'upper left')
    if save_file != "": plt.savefig(save_file,format='png', dpi=300, bbox_inches='tight',transparent=True) 
    plt.show()
    return

def Usage():
    print("Usage:  \n python IFP_contacts_quickView.py IFP_table.pkl \n IFP_table.pkl - pandas frame containing IFPs")


if  len( sys.argv) < 2:
	Usage()
	sys.exit()
finput_file = sys.argv[1]
if not os.path.isfile(finput_file):
        Usage()
        sys.exit()
try:
	IFP_data = pd.read_pickle(finput_file)
except:
        Usage()
        sys.exit()


Plot_IF_trajectory(IFP_data,ifp_type = np.asarray(['AR','HA','HD','IP','IN',"IO","WB"]),head_tail=10,save_file = finput_file[:-3]+".png")
