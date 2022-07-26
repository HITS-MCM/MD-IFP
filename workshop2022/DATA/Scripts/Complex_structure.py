class complex_structure:
    """
    class that contains all informations about the system: a protein/membrane/ligand/water/ion complex
    """
    
       
    def __init__(self,PRJ_DIR = None, PRJ_Name = None):
        """
        Class constructor:
        data will be allocated in dir_complex
        
        """
        self.box = None  #  system box dimention
        self.cyx = None  # set of CYS bridges
        self.PDB = None  # system file in PDB format
        self.ligand_name = [] # a list of residue names that correspond to ligands
    
        self.file_p = None # protein in PDB format
        self.file_w = None # water in PDB format
        self.file_i = None # ions in PDB format
        self.self.file_m = None # membrane in PDB format
    
        self.file_p_aligned = None
        self.file_l_aligned = []
    
        self.ligand_charge = []
        #---------  check/make directory to work in -----------
        if PRJ_DIR:
            if not os.path.exists(PRJ_DIR):
                try:  
                    os.mkdir(PRJ_DIR)
                except OSError:  
                    print ("Creation of the directory %s failed" % PRJ_DIR)
                else:  
                    print ("Data generated will be located at %s " % PRJ_DIR)
        else:
            PRJ_DIR = os.getcwd()
        os.chdir(PRJ_DIR)   
        #-------------------define naming schema----    
        if PRJ_Name:
            pass
        else:
            PRJ_Name = os.path.basename(PRJ_DIR)
        if PRJ_DIR[:-1] != "/": PRJ_DIR = PRJ_DIR+"/"
            
        return

    ##############################################################
    #
    #    In this function a pdb file  prepared by CharmmGUI will be splitted into separate components: protein, membrane, water, ions
    #    Then disulfate bridges will be found and renamed according to the AMBER schema
    #    Finally lipids will be re-constracted by charmmlipid2amber.py
    #    Ligand is not considered at this step    
    #
    ###############################################################
    def charmm_gui(self,charmm_gui_file = "step5_assembly.pdb"):  
        """
        Parameters:
        Results:
        
        """
        from SetUp.setup_membrane_simulations_AMBER import Replace_CYS_CYX , split_charmm_gui
        os.chdir(PRJ_DIR)

        if not os.file.exists(charmm_gui_file):
            print ("File  %s does not exists" % charmm_gui_file)
        else:
            os.chdir(PRJ_DIR)
            #--- split complex generated by Charmm-GUI and find  water box dimention
            self.box = split_charmm_gui(charmm_gui_file,H=False,TER=True)
            #---- find sulfur bridges  and replace CYS by CYX
            self.cyx = Replace_CYS_CYX(charmm_gui_file)
            self.file_p = file_input[:-4]+"_Protein.pdb"
            self.file_i = file_input[:-4]+"_Ions.pdb"
            self.file_m = file_input[:-4]+"_Mem.pdb"
            self.file_w = file_input[:-4]+"_Wat.pdb"
            # --- restructure lipids according to AMBER schema
            line = "charmmlipid2amber.py -i "+self.file_m_name+"  -o  "+self.file_m_name+ "_Amb.pdb"
            proc = subprocess.Popen(line, stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            time.sleep(0.2)
#            proc.terminate()
            try:
                outs, _ = proc.communicate(timeout=0.2)
            except subprocess.TimeoutExpired:
                print('charmmlipid2amber.py did not terminate in time. Check availability of Amber tools')
        return        
             
    ##############################################################
    #
    #    Superimposition of protein (and ligand) with a protein in the system
    #    this will allow us to replace protein in the system by a new protein+ligand
    #
    ###############################################################
    def superimpose_prot_lig(self,protein_pdb,ligand_pdb=[]):
        """
        Parameters:
        Results:
        
        """
        from SetUp.setup_membrane_simulations_AMBER import protein_ligand_align, ligand_clean, protein_ligand_superimpose
        os.chdir(PRJ_DIR)
        self.file_p_aligned = protein_pdb[:-4]+"_aligned.pdb" # new protein PDB filename
        reference_file  = self.file_p_name
        if ligand_pdb:
            self.file_l_aligned = ligand_pdb[:-4]+"aligned.pdb"
            ligand_cleaned = []
            for l in ligand_pdb: 
                ligand_cleaned.append(l[:-4]+"cleaned.pdb")
                check_ligand(l,l[:-4]+"cleaned.pdb")
        protein_ligand_superimpose(reference_file,protein_pdb,self.protein, \
                                   ligand_pdb=ligand_cleaned,new_ligand_pdb = self.file_l_aligned)
        return
        
    ##############################################################
    #
    #    
    #    
    #
    ###############################################################
    def ligand_topology(self,ligand_H_pdb):
        """
        Parameters:
        Results:
        
        """
        from SetUp.setup_membrane_simulations_AMBER import replace_charges, compute_charge
        os.chdir(PRJ_DIR)
        for l in ligand_H_pdb:
            ligand_H_mol2 = l[:-4]+".mol2"
            template_mol2 = l[:-4]+"-template.mol2"
        
            #------- generate charges from the Gamess output----------------
            if(ligand_H_pdb[:-4]+"_energy.dat"):
                self.charge= charge_from_Gamess(ligand_H_pdb[:-4]+"_energy.dat")
            #-----------make a template ligand mol2 file with correct naming according to the Amber schema, but maybe incorrect charges
            if os.file.exists(ligand_pdb):
                if (not self.charge) or (self.charge == 0):  adding = ""
                else:     adding = "-nc "+str(self.charge)+" -m 1 "
                line = os.system(" antechamber -i "+ligand_H_pdb+" -fi pdb -o "+template_mol2+" -fo mol2 -c bcc -s 2 "+adding+ \
                             "; parmchk -i  " + template_mol2 + " -f mol2 -o INH.frcmod") 
                proc = subprocess.Popen(line, stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
                time.sleep(0.2)
            #proc.terminate()
                try:
                    outs, _ = proc.communicate(timeout=0.2)
                except subprocess.TimeoutExpired:
                    print('antechamber did not terminate in time. Check availability of Amber tools')
            #----------- replace charges to the R.E.S.P. values----------------
                if os.path.file(charges.chg):
                    replace_charges(template_mol2,ligand_H_mol2,new_charges_mol2=None,charges=charges.chg)            
            else: print("File  %s was not found" %ligand_H_pdb)
        return