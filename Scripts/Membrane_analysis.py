#!/usr/bin/env python
# coding: utf-8

# # Script for exploring behavior of a membrane-protein system in  MD simulations:
# 
#      1.  membrane/protein/water atom density distribution along z axis (perpendicular to the membrane surface)
#      2.  area per lipid 
#      
#      The class Membrane_properties contains containers for all properties 
#      that can be generated using two main functions of the Membrane_properties class:
#       - Get_info() - computes per-frame system properties
#       - Prep4plot() - computes system properties averaged over all frames 
# 
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
#     trajectory file 
#     pdb file (for example, generated from the first frame)
#     
#     
# ### Packages required:
#     numpy
#     matplotlib
#     MDAnalysis
#     scipy
#     code is written on Python 3.x



get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')
import glob, os
import sys
import numpy as np


from matplotlib import *
from matplotlib import gridspec
import  pylab as plt


import MDAnalysis as mda
from MDAnalysis.lib.formats.libdcd import DCDFile
from MDAnalysis.analysis import contacts,align,rms
from MDAnalysis.analysis.base import AnalysisFromFunction
from MDAnalysis.coordinates.memory import MemoryReader

from scipy.ndimage.filters import gaussian_filter1d
from scipy.interpolate import make_interp_spline, BSpline




#######################################################################
#
#     FUNCTION FOR PUTTING SYSTEM BACK INTO A Periodic Boundary BOX
#
#######################################################################

def pbc(u,Rgr0):
    """
    Parameters:
    u - a frame of the trajectory object (MD universe object)
    Rgr - reference radius of gyration; 
    as a check whether the transformation is correct we compare radius of computed gyration with the reference one
    
    Results:
    Radius of gyration of the protein 
    """
    u_CA = u.select_atoms("name CA")
    sel_p = "protein"
    
    # getting all system elements back to the box; it is important to repeat this twice in the case when protein is splitted into two parts
    u.atoms.translate(-u_CA.center_of_mass()+0.5*u.dimensions[0:3])
    u.atoms.pack_into_box(box=u.dimensions) 
    u.atoms.translate(-u_CA.center_of_mass()+0.5*u.dimensions[0:3])
    u.atoms.pack_into_box(box=u.dimensions) 
    Rgr = u.select_atoms(sel_p).radius_of_gyration()      
    if Rgr > Rgr0*1.1:
#        print("Radius of gyration is too large: %s  of that in the first frame; Try to pack system back into a box once more " %(Rgr/Rgr0)) 
        u.atoms.translate(-u_CA.center_of_mass()+0.5*u.dimensions[0:3])
        u.atoms.pack_into_box(box=u.dimensions) 
        Rgr = u.select_atoms(sel_p).radius_of_gyration()  
#        print("Radius of gyration is now: %s  of the first frame" %(Rgr/Rgr0)) 
    if Rgr > Rgr0*1.1:
        print("failed to pack the system back into a box adius of gyration is too large: %s of that in the first frame" %(Rgr/Rgr0))
    return (Rgr)





class Membrane_properties:
    
    def __init__(self,ref_pdb,sel_ligands = "",interval=(0,-1,1),d=3,dh = 3,sel_m = ["CHL", "PC", "PA", "PE"]):     
        """
        PARAMETERS:
        ref - file with a reference pdb structure
        sel_ligands - ligand residue name (or several names space-separated) to be treated together with tproteion residues
        interval - array containing the first and the last frame, and a stride to be used for analysis
        dz - interval for splitting simulation box (the same in all directions, in Angstrom)
        sel_m - a string defining selection of the residues for membrane analysis. 
            Note, that in AMBER ff PC/PE + PA + OL - are different residues representing one lipid molecules
            particularly, OL and PA are two chains of POPC or POPE lipid that are placed in the same z - interval
            PC or PE are on the top of the lipid and practically do not overlap with PA/OL
            Thus, OL should be omitted to avouid double counting of lipid residues
        d=6 Angstrom over x/y and dh = 3 Angstrom over z 
        Important: 
            - there is no overaging over z, so the best selection for dh is a vdW distance
            - for large x/y steps (dh > 3 A)  the area in x/y plane occupied by ligand/protein will be overestimated
                and atom density will be underestimated because the unit sell is assuming to be completely ocupied by a single atom 
       """
        # input parameters
        self.ref_pdb = ref_pdb
        self.sel_ligands = sel_ligands
        self.interval = interval # defines frames of the trajectory to be processed ( first, last, stirde)
        self.dz = d
        self.dh = dh # approximate vdW distance that defines a unit slab 
        "(not type H ) and ( resname CHL PC PA PE )"
        self.sel_m = "(not type H ) and ( resname "   # for computing lipid number
        self.sel_m_a = "(not type H ) and ( resname "  # for computing all atom number
        for s in sel_m:
            if s != "OL":  self.sel_m = self.sel_m + " " + s 
            self.sel_m_a = self.sel_m_a + " " + s 
        self.sel_m = self.sel_m +" )"
        self.sel_m_a = self.sel_m_a +" )"
        self.margin = max(2,int(4.0/d)+1)  # additional margin used for each side of the box in the case if box size increases along the trajectory
        
        # per-frame data:
            # arrays of numpy 3D matrices: [frame][z, x, y] 
        self.mem_slab = []
        self.prot_slab = []
        self.wat_slab = []
            # arrays of vectors: [frame][z] 
        self.prot_area = []
        self.mem_area = []
        self.wat_area = []
        self.resid_array = []
        self.dens_p = []
        self.dens_w = []
        self.dens_m = []    
        self.m_r_per_area = []
            # arrays of arrays: [frame][z][x]  
        self.resid_array_zx = []
            # vector [frame]
        self.Rgr = []
        
        # data averaged over frames
            # array of numpy matrices [zi][x,y], where zi - includes only z value withing a membrane region (see parameters start_mem,stop.mem, and step.mem)
        self.mem_slab_frame = []
        self.prot_slab_frame = []
            # array of numpy matrices [zi][x], where zi - includes only  membrane region
        self.mem_slab_frame_zx = []
        self.area_x = []
        
        # position of the membrain in the box;  used for plotting
        self.start_mem = None
        self.stop_mem = None
        self.step_mem = None
        return

    ###########################################
    #
    #  Function for computing properties of a
    #  membrane-containing system for a trajectory
    # 
    #############################################

    def Get_info(self,traj):
        """
        PARAMETERS:
        
        traj - trajectory file
        
        Returns:
        
        arrays containing per-frame system analysis:
        1.  arrays of the shape [frame][z, x, y] containing
            mem_slab -  the number of membrane atom
            prot_slab - the number of protein atom 
            wat_slab - the number of water atom 
        2.  arrays of the shape [frame][z] containing
            prot_area  -area occupied by protein atoms
            mem_area   -area occupied by membrane atoms
            wat_area   -area occupied by water atoms
        3.  array of the shape [frame][z] containing
            resid_array - number of lipid resudues (as defined by the parameter sel_m)
        4.  array of the shape [frame][z][x] containing
            resid_array_zx - - number of lipid resudues (as defined by the parameter sel_m)
            for each z,x slab
        Rgr - vector [frame] containing protein radius of gyration for each frame
        """
        ref = mda.Universe(self.ref_pdb)
        Rgr0 = ref.select_atoms("protein").radius_of_gyration() 
        all_atoms = ref.select_atoms("not type H")
    
        # load reference structure
        u = mda.Universe(self.ref_pdb)
        u.load_new(traj)  
        u_length = len(u.trajectory)
        u_size = int(os.path.getsize(traj)/(1024.*1024.))    
        # box parameters
        self.nz = int(u.dimensions[2]/self.dh)+1+self.margin
        self.nx = int(u.dimensions[0]/self.dz)+self.margin
        self.ny = int(u.dimensions[1]/self.dz)+self.margin
        print("DIM (from traj): ",u.dimensions)
        print("DIM (x/y/z): ",self.nx,self.ny,self.nz)
        
        u_mem = u  # can be replaced by a procedure of loading trajectory in RAM
        u_mem_length = u_length
        start_analysis = nx = self.interval[0] 
        stop_analysis = self.interval[1] 
        step_analysis = self.interval[2] 
        if(stop_analysis < 0): stop_analysis = u_mem_length
        system_reduced = u.select_atoms(" protein ")
        frames = int((stop_analysis-start_analysis)/step_analysis)
        nx = self.nx
        ny = self.ny
        nz = self.nz
        print("number of frames= %s; file size %s M" %(u_length,u_size))
        print("will be analyzed  %s frames" %(frames))

        wat_slab = []
        #sel_m = "(not type H Cl Na NA CL Cl- Na+) and (not protein) and (not resname WAT "+sel_ligands+" 2CU)  "
        #sel_m = "(not type H ) and ( resname CHL PC PA PE )  "  # PA + OL + PC makes POPC; PA and OL are two tails of POPC
        sel_m = self.sel_m
        sel_m_a = self.sel_m_a
        sel_w = "(not type H ) and(resname WAT)"
        if len(self.sel_ligands) > 1:
            sel_p = "(not type H) and  (protein  or  (resname "+self.sel_ligands+"))"
        else:
            sel_p = "(not type H) and  protein" 
    
        self.Rgr = []
        for i0,i in enumerate(range(start_analysis,stop_analysis,step_analysis)):
            if (i0%10 == 0):
                print("frames analized: ",i0," current frame: ",i)
            u_mem.trajectory[i]
            self.Rgr.append(pbc(u_mem,Rgr0))
            u_mem_sel_m = u_mem.select_atoms(sel_m)
            u_mem_sel_m_a = u_mem.select_atoms(sel_m_a)
            u_mem_sel_p = u_mem.select_atoms(sel_p)            
            u_mem_sel_w = u_mem.select_atoms(sel_w)                    
            
            mem_slab0 = np.zeros((nz,nx,ny),dtype = int)
            resid_list = []
            resid_list_zx = []
            for t in range(0,nz): 
                resid_list.append([])
                resid_list_zx.append([])
                for t in range(0,nx): resid_list_zx[-1].append([])
            if(sel_m_a != sel_m):
                for at,t in zip(u_mem_sel_m.positions,u_mem_sel_m):
                    n = t.resid
                    ix = (int)(at[0]/self.dz)
                    iy =  (int)(at[1]/self.dz)
                    iz = (int)(at[2]/self.dh)
                    try:
                        mem_slab0[iz,ix,iy] += 1
                    except:
                        print("Membrane is out of the box ",at," ",nx,ny,nz," ",ix,iy,iz)
                for at,t in zip(u_mem_sel_m_a.positions,u_mem_sel_m_a):
                    n = t.resid
                    ix = (int)(at[0]/self.dz)
                    iy =  (int)(at[1]/self.dz)
                    iz = (int)(at[2]/self.dh)
                    try:
                        if n not in resid_list[iz]: 
                            resid_list[iz].append(n)
                            resid_list_zx[iz][ix].append(n)
                    except:
                        pass
            else:
                for at,t in zip(u_mem_sel_m.positions,u_mem_sel_m):
                    n = t.resid
                    ix = (int)(at[0]/self.dz)
                    iy =  (int)(at[1]/self.dz)
                    iz = (int)(at[2]/self.dh)
                    try:
                        mem_slab0[iz,ix,iy] += 1
                        if n not in resid_list[iz]: 
                            resid_list[iz].append(n)
                            resid_list_zx[iz][ix].append(n)
                    except:
                        print("Membrane is out of the box ",at," ",nx,ny,nz," ",ix,iy,iz)
               
            prot_slab1 = np.zeros((nz,nx,ny),dtype = int)
            for at,t in zip(u_mem_sel_p.positions,u_mem_sel_p):
                    ix = (int)(at[0]/self.dz)
                    iy =  (int)(at[1]/self.dz)
                    iz = (int)(at[2]/self.dh)
                    try:
                        prot_slab1[iz,ix,iy] += 1
                    except:
                        print("Protein is out of the box ",at," ",nx,ny,nz," ",ix,iy,iz)
                        
            wat_slab2 = np.zeros((nz,nx,ny),dtype = int)
            for at,t in zip(u_mem_sel_w.positions,u_mem_sel_w):
                    ix = (int)(at[0]/self.dz)
                    iy = (int)(at[1]/self.dz)
                    iz = (int)(at[2]/self.dh)
                    try:
                        wat_slab2[iz,ix,iy] = 1
                    except:
                        if iz < nz:  # wather at box boundary is not important, just skip it
                            print("Water is out of the box ",at," ",nx,ny,nz," ",ix,iy,iz)
                        
            mem_area1 = np.zeros((nz),dtype = int)
            prot_area1 = np.zeros((nz),dtype = int)
            resid_array0 = np.zeros((nz),dtype = int)
            resid_array_zx0 = np.zeros((nz,nx),dtype = int)
            wat_area1 = np.zeros((nz),dtype = int)
            #---- ToDo - computing of the area can be done more accurately by taking into account just the area occupied by a single atom
            for iz in range(0,nz): 
                """
                prot_area1[iz] =  np.sum(prot_slab1[iz]) *4.0
                mem_area1[iz] =  np.sum(mem_slab0[iz]) *4.0
                wat_area1[iz] =  np.sum(wat_slab2[iz]) *4.0
                """
                prot_area1[iz] = np.count_nonzero(prot_slab1[iz])*self.dz*self.dz
                mem_area1[iz] = np.count_nonzero(mem_slab0[iz])*self.dz*self.dz
                wat_area1[iz] = np.count_nonzero(wat_slab2[iz])*self.dz*self.dz
                resid_array0[iz] = len(resid_list[iz])
                for ix in range(0,nx):
                    resid_array_zx0[iz][ix] = len(resid_list_zx[iz][ix])
            self.wat_area.append(wat_area1)
            self.mem_area.append(mem_area1)
            self.prot_area.append(prot_area1)
            self.mem_slab.append(mem_slab0)
            self.wat_slab.append(wat_slab2)
            self.prot_slab.append(prot_slab1)
            self.resid_array.append(resid_array0)
            self.resid_array_zx.append(resid_array_zx0)
        return 
    
    
   
    ##################################################
    #
    # Averagion membrane properties over all frames and
    # Preparation of data for plotting
    #
    ##################################################
    def Prep4plot(self):
        """
        Parameters:
        
        Results:
            dens_p - density of  protein (and ligand) atoms
            dens_m - density of  membrane atoms
            m_r_per_area - membrane residues per squered Angstrom
            
            mem_slab_frame - [z][x,y] x/y distribution of membrane atoms for a set of z slabs
            prot_slab_frame - [z][x,y] x/y distribution of protein atoms for a set of z slabs
            mem_slab_frame_zx - [z][x] number of membrane  residues for particular z and x (summed up over y)
            area_x - [z][x]  area occupied by membrane atoms for particular z and x (summed up over y)
        
        """
  
        area_xy = self.nx*self.ny*self.dz*self.dz

        # first we will estimate the position of the membrane center
        self.start_mem = np.argwhere(self.resid_array[0] > 0)[0][0]+1
        self.stop_mem = np.argwhere(self.resid_array[0] > 0)[-1][0]-1
        self.step_mem = 2

        area_p = self.prot_area
        frames = len(self.resid_array)
    
        for frame in range(0,frames):
            # ---- density of protein atoms  as function of z
            number_p = np.sum(np.sum(self.prot_slab[frame], axis=2),axis=1)
            dens_p0 = []
            for n,(a,d) in enumerate(zip(area_p[frame],number_p)):
                if(a > 0 ): dens_p0.append(d/a)
                else: dens_p0.append(0)           
            self.dens_p.append(np.asarray(dens_p0))
    
            area_m = (area_xy-self.prot_area[frame]) #-wat_area[frame])
            number_m = np.sum(np.sum(self.mem_slab[frame], axis=2),axis=1)
            dens_m0 = []
            dens_r0 = []
    
            # ---- density of membrane atoms and residues as function of z   
            for n,(a,d,r) in enumerate(zip(area_m,number_m,self.resid_array[frame])): # loop over z
                # density is non-zero only if area is non-zero, the number of lipids is more than 10 and the number of atoms is more than 75
                if(a > 0 ):  #and r > 75 d > 10
                    dens_m0.append(d/a)
                    if(d/a > 0.025):  dens_r0.append(a/r)
                    else:  dens_r0.append(0)
                else: 
                    dens_m0.append(0)
                    dens_r0.append(0)


            self.dens_m.append(np.asarray(dens_m0))
            self.m_r_per_area.append(np.asarray(dens_r0))
        
            # ---- density of membrane residues as function of z and x
            for i,z in enumerate(range(self.start_mem ,self.stop_mem,self.step_mem)): # loop over z and average over frames
                if(frame == 0): 
                    self.mem_slab_frame.append(self.mem_slab[frame][z])
                    self.prot_slab_frame.append(self.prot_slab[frame][z])  
                    self.area_x.append(self.dz*self.dz*(self.nx-np.count_nonzero(self.prot_slab[frame][z],axis=1)))
#                    self.area_x.append(self.dz*self.dz*self.nx-4*np.sum(self.prot_slab[frame][z],axis=1))
                    self.mem_slab_frame_zx.append(self.resid_array_zx[frame][z])  # dencity of lipids
                else: 
                    self.mem_slab_frame[i] = np.add(self.mem_slab_frame[i],self.mem_slab[frame][z])
                    self.prot_slab_frame[i] = np.add(self.prot_slab_frame[i],self.prot_slab[frame][z])
#                    self.area_x[i]  = np.add(self.area_x[i],self.dz*self.dz*self.nx-4*np.sum(self.prot_slab[frame][z],axis=1))
                    self.area_x[i]  = np.add(self.area_x[i],self.dz*self.dz*(self.nx-np.count_nonzero(self.prot_slab[frame][z],axis=1)))
                    self.mem_slab_frame_zx[i] = np.add(self.mem_slab_frame_zx[i] ,self.resid_array_zx[frame][z])
        return

    ###########################################
    #
    # check if vectors for xz is in agreement for vectors for z;
    #
    ###########################################
    def Check(self):
        plt.plot(self.resid_array[0], label="lipids")
        plt.plot(self.prot_area[0], label="prot. area")

        plt.plot(np.sum(self.resid_array_zx[0],axis=1), alpha = 0.2,lw = 10,  label="lipids from xy")
        plt.plot(np.sum(np.count_nonzero(self.prot_slab[0],axis=1),axis=1), alpha = 0.2,lw =10, label="protein atoms")
        plt.legend(framealpha = 0.0,edgecolor ='None',loc='best')
        plt.ylabel('arb.un.', fontsize=14)
        plt.xlabel('z', fontsize=14)
        plt.show()
        return

    ###########################################
    #
    # Plot of 
    #   1. protein/membrane/water atom density as a function of z (axis perpendicular to a membrane surfce)
    #   2. area per lipid  as a function of z
    #
    ###########################################
    def Plot_mem_prot_wat_dens(self):
        """
        PARAMETERS:
        dens_p,dens_m - density of the protein and membrane atoms, respectively, as a function of z
        m_r_per_area - membrane residues per squered Angstrom  as a function of z
        prot_area, mem_area,wat_area - area occupied by atoms of protein, membrane, and water, respectively, as a function of z

        """
    
        dens_p = self.dens_p
        dens_m = self.dens_m
        m_r_per_area = self.m_r_per_area
        prot_area = self.prot_area
        mem_area = self.mem_area
        wat_area = self.wat_area
    
        X = self.dh*np.asarray(range(0,len(dens_p[0])))            
        fig = plt.figure(figsize=(12, 6),dpi=150)
        gs = gridspec.GridSpec(2, 2,hspace=0.5) #height_ratios=[2,2,1]) #,width_ratios=[2,2,1,1])
        ax2 = plt.subplot(gs[0])
        plt.errorbar(x=X,y=np.mean(np.asarray(dens_p),axis=0),  yerr= np.std(np.asarray(dens_p),axis=0), color = "gray" , fmt='o--', markersize=1)
        plt.scatter(x=X,y=np.mean(np.asarray(dens_p),axis=0),color = 'red',alpha=0.5,s=50, label="protein")

        plt.errorbar(x=X,y=np.mean(np.asarray(dens_m),axis=0), yerr= np.std(np.asarray(dens_m),axis=0), color = "gray" , fmt='o--', markersize=1 )
        plt.scatter(x=X,y=np.mean(np.asarray(dens_m),axis=0),color = 'green',alpha=0.5,s=50, label="membr.")

        ax2.legend(framealpha = 0.0,edgecolor ='None',loc='best')
        ax2.set_ylabel('[atoms/A^2]', fontsize=14)
        ax2.set_xlabel('z-distance', fontsize=14)
        ax2.set_title('Density ')
        ax2.grid(color='gray', linestyle='-', linewidth=0.2)

        ax4 = plt.subplot(gs[1])

        plt.errorbar(x=X,y=np.mean(np.asarray(m_r_per_area),axis=0),yerr= np.std(np.asarray(m_r_per_area),axis=0), color = "gray" , fmt='o', markersize=1)
        plt.scatter(x=X,y=np.mean(np.asarray(m_r_per_area),axis=0),color = 'green',alpha=0.5,s=50, label="membr.")
        ysmoothed = gaussian_filter1d(np.mean(np.asarray(m_r_per_area),axis=0), sigma=1)
        plt.plot(X,ysmoothed,color = "green" , lw = 1)
        ax4.set_ylim(0,max(np.mean(np.asarray(m_r_per_area),axis=0)))
        ax4.legend(framealpha = 0.0,edgecolor ='None',loc='best')
        ax4.set_ylabel('area [A^2]', fontsize=14)
        ax4.set_xlabel('z-distance', fontsize=14)
        ax4.set_title('Area per lipid')
        ax4.grid(color='gray', linestyle='-', linewidth=0.2)

        ax3 = plt.subplot(gs[2])
        plt.errorbar(x=X,y=np.mean(np.asarray(prot_area),axis=0),             yerr= np.std(np.asarray(prot_area),axis=0), color = "gray" , fmt='o--', markersize=1)
        plt.scatter(x=X,y=np.mean(np.asarray(prot_area),axis=0),color = 'red',alpha=0.5,s=50, label="protein")

        plt.errorbar(x=X,y=np.mean(np.asarray(mem_area),axis=0),             yerr= np.std(np.asarray(mem_area),axis=0), color = "gray" , fmt='o--', markersize=1 )
        plt.scatter(x=X,y=np.mean(np.asarray(mem_area),axis=0),color = 'green',alpha=0.5,s=50, label="membr. ")

        plt.errorbar(x=X,y=np.mean(np.asarray(wat_area),axis=0),             yerr= np.std(np.asarray(wat_area),axis=0), color = "gray" , fmt='o--', markersize=1 )
        plt.scatter(x=X,y=np.mean(np.asarray(wat_area),axis=0),color = 'blue',alpha=0.5,s=50, label="water ")

        ax3.legend(framealpha = 0.0,edgecolor ='None',loc='best')
        ax3.set_ylabel(r' $Angstrom^2 $ ', fontsize=14)
        ax3.set_xlabel('z-distance [A]', fontsize=14)
        ax3.set_title('Area in the xy plane occupied by atoms of each part of the system')
        ax3.grid(color ='gray', linestyle='-', linewidth=0.2)
        
        ax5 = plt.subplot(gs[3])
        plt.scatter(x = self.interval[2]*np.asarray(range(0,len(self.Rgr))),y=self.Rgr,color = 'blue',alpha=0.5,s=50)
        ax5.set_ylabel('Rad. of gyration [A]', fontsize=14)
        ax5.set_xlabel('frame', fontsize=14)
        ax5.set_title('Rad. of Gyration (protein) ')
        ax5.grid(color='gray', linestyle='-', linewidth=0.2)
        ax5.set_ylim(0.5*max(self.Rgr),1.2*max(self.Rgr))
        plt.show()  
        return


    ###########################################
    #  
    # Plot of 
    #    1. atom distribution for a protein and a membrane at several z distances
    #    2. area per lipid at several z and x distances (averaged over y)
    #
    ###########################################

    def Plot_mem_z(self):
        """
        PARAMETERS:
        prot_slab_frame,mem_slab_frame - (z,x,y) matrix of atom densities  for protein and membrane, respectively
        mem_slab_frame_zx -(z,x) membrane residue dencity in the x,z slab
        area_x - (x,z) matrix of protein occupied area 
   
        """
        prot_slab_frame = self.prot_slab_frame
        mem_slab_frame = self.mem_slab_frame
        mem_slab_frame_zx = self.mem_slab_frame_zx
        area_x = self.area_x

        start_mem = self.start_mem
        stop_mem = self.stop_mem
        step_mem = self.step_mem
    
        pl = 0
        dens_m_xy = []
        fig = plt.figure(figsize=(16, 7))
        plots = len(prot_slab_frame)
        gs = gridspec.GridSpec(4, plots,hspace=0.1,wspace = 0.5) #,height_ratios=[1,1],width_ratios=[2,2,1,1])
        middle = int(0.5*len(mem_slab_frame_zx))
        zmax = 1.5*np.max(np.divide(area_x[middle][self.margin:-self.margin],mem_slab_frame_zx[middle][self.margin:-self.margin]))
    
        for i,z in enumerate(range(start_mem ,stop_mem,step_mem)):
            ax1 = plt.subplot(gs[pl])    
            ax1.imshow(prot_slab_frame[i], interpolation='hamming', cmap="Reds")
            ax1.set_title('z='+str(self.dh*z))
            ax1.set_ylabel('protein', fontsize=16)
            plt.yticks([])
            plt.xticks([])
   
            ax3 = plt.subplot(gs[pl+plots])
            ax3.imshow(mem_slab_frame[i], interpolation='hamming', cmap="Greens")
            ax3.set_ylabel('membrane', fontsize=16)
            plt.yticks([])
            plt.xticks([])

            ax4 = plt.subplot(gs[pl+2*plots])
            Y = np.divide(area_x[i][self.margin:-self.margin],mem_slab_frame_zx[i][self.margin:-self.margin])
            plt.scatter(x=self.dz*np.asarray(range(0,Y.shape[0])),y=Y,color = 'green',alpha=0.5,s=30)
            ysmoothed = gaussian_filter1d(Y, sigma=1)
            plt.plot(self.dz*np.asarray(range(0,Y.shape[0])),ysmoothed,color = "green" , lw = 3,alpha=0.5)
            ax4.set_ylim(0,zmax)
            ax4.grid(color='gray', linestyle='-', linewidth=0.2)
            if(i == 0): 
                ax4.set_ylabel('area/lipid [A^2]', fontsize=16)
                ax4.set_xlabel('x-distance', fontsize=16)

            pl += 1
        plt.show()  
        return








