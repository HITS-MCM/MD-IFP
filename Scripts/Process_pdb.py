import re
import wget
import os
import subprocess
import MDAnalysis as mda
import warnings
import pandas as pd
from Bio.PDB import *
from Bio import BiopythonWarning
from Scripts.IFP_generation import IFP, table_combine
from rdkit import Chem, RDConfig


def extract_conformation(line, conformation, ligand):
    """Helper function to extract information about conformation from pdb file.

    Beware of bugs. pdb format is very inconsistent and columns might merge together.


    Args:
        line (str): Line from pdb file that will get analyzed.
        conformation (str): Conformation that will be found.
        ligand:

    Returns:
        1 if conformation information is matched.
        0 if conformation information is not matched.
    """
    residues = standard_aa_names + [ligand]
    split_line = line.split()
    if split_line[0][:6].strip() in ["ATOM", "HETATM"]:
        if split_line[3][-3:] in residues:
            try:
                if split_line[3][-4] in conformation:
                    return 1
                else:
                    return 0
            except IndexError:
                return 0
        elif split_line[2][-3:] in residues:
            try:
                if split_line[2][-4] in conformation:
                    return 1
                else:
                    return 0
            except IndexError:
                return 0

        elif split_line[1][-3:] in residues:
            try:
                if split_line[1][-4] in conformation:
                    return 1
                else:
                    return 0
            except IndexError:
                return 0

    else:
        return 0


class PdbIDAnalysis(object):
    """Prepares files from PDB-ID using Chimera and automatically performs IFP analysis.

    Crystal structures are downloaded from rcsb.org, split into the given protein chain, ligand and water molecules.
    Using Chimera hydrogens are added and mol2 of the ligand is generated. The three pdb files are merged. IFP Analysis
    is then performed.
    If the file contains multiple conformations IFPs will not be found correctly. The conformations need to be filtered
    out manually using the split_pdb method.

    Class attributes:
        reading_failed (list): PDB-IDs for which rdkit was unable to properly read the ligand.mol2 file.
        features_failed (list): PDB-IDs for which rdkit was unable to create features.
    """
    features_failed = []
    reading_failed = []

    def __init__(self, pdb_id, ligand_id, data_path='../Data/', verbose=False):
        """Inits PdbIDAnalysis with given parameters.

        Args:
            pdb_id (str): PDB-ID for which the analysis will be performed.
            ligand_id (str): ID of the ligand which is supposed to be analyzed.
            data_path (str): Path indicating where generated data is supposed to be stored. Defaults to parent directory.
            verbose (bool): If true prints additional information during the processing. Defaults to false.
        """
        if not os.path.isdir(data_path):
            os.mkdir(data_path)
        self.data_path = data_path + pdb_id + "/"
        self.verbose = verbose
        self.pdb_id = pdb_id
        self.ligand_id = ligand_id
        if not os.path.isdir(self.data_path):
            os.mkdir(self.data_path)

    # Downloads the given pdb dataset from the RCSB. It checks for correct formatting and folder structure
    def download_pdb(self):
        """ Downloads pdb file from rcsb.org.

        For downloading uses the URL supplied by rcsb.org for uncompressed pdb files.
        If given data path is not existing it is created. If the file already exists, download is not performed.

        Raises:
            ValueError: The supplied PDB-Id does not conform to the format standard.
        """

        # Checks that pdb_id starts with a digit and is not longer than 4 characters.
        if (re.match(r"[1-9]\w\w\w", self.pdb_id) is None) or len(self.pdb_id) != 4:
            raise ValueError('Input data does not conform to pdb standard')

        file_name = self.data_path + self.pdb_id + '.pdb'

        if os.path.isfile(file_name) is True:
            if self.verbose is True:
                print('File {} has already been downloaded.'.format(self.pdb_id))

        else:
            if self.verbose is True:
                print('Downloading {}.pdb from RCSB to {}'.format(self.pdb_id, file_name))
            wget.download('https://files.rcsb.org/download/{0}.pdb'.format(self.pdb_id), out=file_name)
        return

    def split_pdb(self, pdb_chain, unwanted_conformation=None, wanted_conformation=None, residue_id = None):
        """Extracts specified chain, ligand and water molecules from pdb file.

        Sub-structures are extracted using the pdb-handling capacities of the Biopython package. A new folder is created
        for the split structures.
        The with ligand_conformation specified conformation will be filtered out.

        Args:
            pdb_chain (str): Chain identifier for the chain that is supposed to be extracted.
            unwanted_conformation (str):  Conformation to be filtered out. Not mentioned conformation remains in pdb.
            wanted_conformation (str): Conformation that is wanted. Used to reformat lines in case atom description and
                                       residue merge.
            residue_id (int): Optional wanted residue number of the ligand. Usefull if two ligands interact with the
                              same chain.

        Returns:
            None, if pdbID_split directory already exists.
        """
        if unwanted_conformation is not None and not os.path.isfile(
                self.data_path + "unfiltered_" + self.pdb_id + ".pdb"):
            os.rename(self.data_path+self.pdb_id+".pdb", self.data_path+"unfiltered_" + self.pdb_id + ".pdb")
            with open(self.data_path + "unfiltered_" + self.pdb_id + ".pdb", "r+") as file:
                lines = file.readlines()
                with open(self.data_path + self.pdb_id + ".pdb", "w+") as filtered:
                    for line in lines:
                        if not extract_conformation(line, unwanted_conformation, self.ligand_id):
                            if extract_conformation(line, wanted_conformation, self.ligand_id):
                                fixed_line = ""
                                for residue in standard_aa_names + [self.ligand_id]:
                                    index = line.find(residue)
                                    if index == -1:
                                        continue
                                    else:
                                        fixed_line = fixed_line + line[:index-1] + " "
                                        fixed_line = fixed_line + line[index:]
                                filtered.write(fixed_line)
                            else:
                                filtered.write(line)

        structure = PDBParser().get_structure(self.pdb_id, self.data_path + self.pdb_id + '.pdb')
        ligand_id = self.ligand_id
        split_path = self.data_path + self.pdb_id + "_split/"
        if not os.path.isdir(split_path):
            os.mkdir(split_path)

        else:
            if self.verbose is True:
                print("File {0} seems to have already been split.".format(self.pdb_id))
            return

        class LigandSelect(Select):
            """Selects only atoms belonging to the given ligand ID.

            accept_residue is overloaded to only return 1 if the residue name is equal to the ligand ID.
            If residue_id is specified 1 will only be returned if the residue_id is matched.
            accept_chain is overloaded to only return 1 if the chain is equal to the supplied chain identifier.

            """
            def accept_residue(self, residue):
                if residue.get_resname().strip() == ligand_id:
                    if residue_id is not None:
                        if residue.get_id()[1] == residue_id:
                            return 1
                        else:
                            return 0
                    else:
                        return 1
                else:
                    return 0

            def accept_chain(self, chain):
                if chain.get_id() == pdb_chain:
                    return 1
                else:
                    return 0

        class WaterSelect(Select):
            """Selects only atoms belonging to water molecules.

            accept_residue is overloaded to only return 1 if the residue name is HOH.
            accept_chain is overloaded to only return 1 if the chain is equal to the supplied chain identifier.

            """
            def accept_residue(self, residue):
                if residue.get_resname() == "HOH":
                    return 1
                else:
                    return 0

            def accept_chain(self, chain):
                if chain.get_id() == pdb_chain:
                    return 1
                else:
                    return 0

        class ChainSelector(Dice.ChainSelector):
            """Selects only atoms belonging to the given protein chain.

            accept_residue is overloaded to only return 1 if the residue is not flagged as a Heteroatom. Compared to the
            parent method length limitation has been removed.

            accept_atom is overloaded to return 1 for every atom type. Parent method only returns 1 for non-hydrogen
            atoms.

            """
            def __init__(self, chain_id, model_id=0):
                self.chain_id = chain_id
                self.model_id = model_id

            def accept_residue(self, residue):
                hetatm_flag, resseq, icode = residue.get_id()
                if hetatm_flag != " ":
                    return 0
                if icode != " ":
                    warnings.warn("WARNING: Icode %s at position %s" % (icode, resseq), BiopythonWarning)
                return 1

            def accept_atom(self, atom):
                return 1

        io = PDBIO()
        io.set_structure(structure)

        if self.verbose is True:
            print("Extracting protein chain {0}, ligand {2} and water. Saving to {1}".format(pdb_chain, split_path,
                                                                                             self.pdb_id))
        io.save(split_path + "protein" + '.pdb', ChainSelector(pdb_chain))
        io.save(split_path + self.ligand_id + '.pdb', LigandSelect())
        io.save(split_path + 'water' + '.pdb', WaterSelect())

        with open(self.data_path + self.pdb_id + "_split/merged.pdb", "w+") as outfile:
            for file in ["protein", self.ligand_id, "water"]:
                with open(self.data_path + self.pdb_id + "_split/" + file + ".pdb", "r") as infile:
                    # Only Atom/Hetatom lines are written
                    lines = [line for line in infile if line.split()[0] not in ["CONECT", "MASTER", "END", "TER"]]
                    outfile.writelines(lines)

    def add_hydrogen(self, script="add_hydrogen_chimera.py",
                     chimera_path="C:/Program Files/Chimera 1.14/bin/chimera.exe"):
        """chimera_hydrogen_mol2.py is executed in Chimera and pdb files are merged.

        Through the shell Chimera is opened without gui and chimera_hydrogen_mol2.py is run with the split pdb files as
        input. The three edited pdb files are then merged into one called merged.pdb

        Args:
            script (str): Script that chimera is supposed to run.
            chimera_path (str): Path to the chimera executable.
        """
        if os.path.isdir(self.data_path + self.pdb_id + "_split/"):
            return
        subprocess.run("\"{0}\" --nogui --script \"{1}/{2} {3}\"".format(
            chimera_path, os.getcwd(), script, self.data_path), shell=True)

    def create_ligand_mol2(self, script="create_mol2_chimera.py",
                           chimera_path="C:/Program Files/Chimera 1.14/bin/chimera.exe"):
        if os.path.isfile(self.data_path + self.pdb_id + "_split/" + self.ligand_id + ".mol2"):
            return
        subprocess.run("\"{0}\" --nogui --script \"{1}/{2} {3} {4}\"".format(
            chimera_path, os.getcwd(), script, self.data_path+self.pdb_id+"_split", self.ligand_id), shell=True)

    def calculate_ifp(self, get_properties=False, wanted_conformation=None):
        """ IFPs are calculated.

        Properties of the ligand are determined using ligand mol2 and pdb files. merged.pdb is loaded into an MDAnalysis
        Universe. Properties and the Universe are used to compute the IFPs, which are combined and transformed into a
        dataframe.

        Args:
            get_properties (bool): If true additionally returns the ligand properties.
            wanted_conformation (str): Name of the wanted conformation if multiple conformations in one file.

        Returns:
            Pandas Dataframe containing the IFPs.

        """
        ligand_pdb = self.data_path + self.pdb_id + "_split/" + self.ligand_id + ".pdb"
        ligand_mol2 = self.data_path + self.pdb_id + "_split/" + self.ligand_id + ".mol2"
        pdb_merged = self.data_path + self.pdb_id + "_split/merged.pdb"

        with open(ligand_pdb, 'r') as fasta_file:
            if wanted_conformation is None:
                list_labels = [line.split()[2] for line in fasta_file.readlines() if line.split()[0] in ["ATOM", "HETATM"]]
            elif wanted_conformation:
                list_labels = []
                for line in fasta_file.readlines():
                    if line.split()[0] in ["ATOM", "HETATM"]:
                        if len(line.split()) == 11:
                            list_labels.append(line.split()[2][:-4])
                        else:
                            list_labels.append(line.split()[2])

        try:
            mol = Chem.rdmolfiles.MolFromMol2File(ligand_mol2, removeHs=False)

        # If an error occurs during reading of ligand.mol2, the pdb ID is appended to reading_failed and skipped.
        except OSError:
            warnings.warn("mol2 file of the ligand {0} could not be read.".format(self.ligand_id))
            PdbIDAnalysis.reading_failed.append(self.pdb_id)
            return

        factory_name = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
        factory = Chem.ChemicalFeatures.BuildFeatureFactory(factory_name)

        try:
            features = factory.GetFeaturesForMol(mol)

        # If an error occurs during feature generation of the ligand, the pdb ID is appended to features_failed and
        # skipped.
        except:
            warnings.warn("Features for the ligand {0} could not computed".format(self.ligand_id))
            PdbIDAnalysis.features_failed.append(self.pdb_id)
            return

        properties = {feature.GetFamily(): [] for feature in features}
        for feature in features:
            current_atoms = [list_labels[atom_id] for atom_id in list(feature.GetAtomIds())]
            for atom in current_atoms:
                properties[feature.GetFamily()].append(atom)

        u = mda.Universe(pdb_merged)

        df_prop, df_hb, df_wb = IFP(u, self.ligand_id, properties, RE=False)
        df_prop_complete = table_combine(df_hb, df_wb, df_prop, self.ligand_id)
        df_prop_complete.rename(columns={"time": self.ligand_id}, inplace=True)

        if get_properties is True:
            return df_prop_complete, properties

        return df_prop_complete

    def run(self, chain, chimera_path="C:/Program Files/Chimera 1.14/bin/chimera.exe", get_properties=False,
            residue_id=None):
        """ Wraps all methods of PdbIDAnalysis into a single method for convenient use. For Args see each method itself.

        """
        self.download_pdb()
        self.add_hydrogen(chimera_path=chimera_path)
        self.split_pdb(chain, residue_id=residue_id)
        try:
            self.create_ligand_mol2(chimera_path=chimera_path)
        except:
            PdbIDAnalysis.features_failed.append(self.pdb_id)
        return self.calculate_ifp(get_properties=get_properties)


def merge_hd_with_ha(dataframe):
    """Merges hydrogen donor column with hydrogen acceptor column.

    Since in the validation set there was no differentiation in hydrogen bonds between hydrogen donor and
    hydrogen acceptor and our IFP genereation does differentiate, both interactions are combined for easier analysis.

    Args:
        dataframe (pd.DataFrame): Interaction Dataframe with HD and HA column that will be merged.

    Returns:
        Dataframe with HD and HA columns merged into HD column.
    """
    for protein in dataframe.index:
        buffer = dataframe.loc[protein, "HA"]
        for residue in buffer:
            dataframe.loc[protein, "HD"].append(residue)
    return dataframe.drop(columns="HA")


class StatisticalAnalysis(object):
    """Compares results with baseline.

    Baseline is taken from Hajiebrahimi et al., 2017.

    """
    def __init__(self, data, interactions=("AR", "HD", "HA", "HL", "IN", "IP", "WB")):
        """Inits StatisticalAnalysis with entered data.

        Args:
            data (str): Path to generated IFPs saved in pickle format.
            interactions (tuple): Interactions for which comparison should be done.
        """
        self.data = pd.read_pickle(data)
        self.possible_interactions = list(interactions)

    def get_results(self):
        """Gets generated IFPs into analyzable format.

        Returns:
            pandas Dataframe object.
            Columns: Elements of possible_interactions attribute.
            Index: PDB IDs of structures.
            Values: Residue Identifiers that build the interaction mentioned in the column.

            First column is named "ligand" and contains the ligand ID corresponding to the PDB ID.

        """

        def extract_interactions(pdb_id, interaction_list):
            """Helper function that extracts residue in IFP and assigns an interaction.

            Splits entrys from format interaction_residue into interaction and residue. Residues are collected
            according to their interaction.

            Args:
                pdb_id (str): PDB ID. Used to index and get entrys from the input DataFrame.
                interaction_list (list): Possible interactions.

            Returns:
                A Dictionary mapping interaction type to residues that interact with the ligand in that type.
                Example:
                    {"Aromatic": ("Phe1", "Tyr2", "His3"),
                    "Ionic": ("Arg4", "His5", "Lys6")}
            """
            interaction_dict = {element: [] for element in interaction_list}
            for entry in self.data.loc[pdb_id]:
                if entry not in [None, "time", "WAT"]:
                    entry_split = entry.split("_")
                    try:
                        interaction_dict[entry_split[0]].append(entry_split[1])
                    # If interactions is not supposed to be found skip it.
                    except KeyError:
                        continue
                else:
                    continue
            return interaction_dict

        results = []
        for index in self.data.index:
            results.append(extract_interactions(index, self.possible_interactions))
        results = pd.DataFrame(results, self.data.index)
        results.insert(0, "ligand", self.data[0])
        return results

    def compare(self, baseline, results):
        """False positives and false negatives are computed.

        Args:
            baseline (pd.DataFrame): Used as ground truth.
            results (pd.DataFrame): Used as observation.

        Returns:
            Dataframe containing residues appearing in observation but not in the ground truth (false postivies) broken
            down after interaction and PDB ID.

            Dataframe containing residues appearing in ground truth but not in the observation (false negatives) broken
            down after interaction and PDB ID.
        """

        baseline = merge_hd_with_ha(baseline)
        results = merge_hd_with_ha(results)

        false_positives = pd.DataFrame([], index=self.data.index, columns=self.possible_interactions)
        false_negatives = pd.DataFrame([], index=self.data.index, columns=self.possible_interactions)

        for protein in results.index:
            for interaction in self.possible_interactions:
                if interaction == "HA":
                    continue

                fal_pos = [residue for residue in results.loc[protein, interaction]
                           if residue not in baseline.loc[protein, interaction]]
                false_positives.loc[protein, interaction] = fal_pos

                fal_neg = [residue for residue in baseline.loc[protein, interaction]
                           if residue not in results.loc[protein, interaction]]
                false_negatives.loc[protein, interaction] = fal_neg

        return false_positives.drop(columns="HA"), false_negatives.drop(columns="HA")

    def get_baseline(self, input_file):
        """ Placeholder method that is used since baseline is badly formatted and can't easily be used.
        """

        result = {}

        with open(input_file, "r") as input_fil:
            data = input_fil.read()

            for lines in data.split("\n"):
                line = lines.split("\t")

                if len(line[0]) != 4:
                    continue

                entry = [[residue.strip() for residue in interaction.split(",")] for interaction in line[2:]]
                result[line[0]] = entry

        return pd.DataFrame(result, index=self.possible_interactions).T


if __name__ == '__main__':
    hajiebrahimi_pdb = ["1bju", "1bma", "1eve", "1hvy", "1hwi", "1ia1", "1ig3", "1jd0", "1jla", "1l2s",
                        "1l7f", "1n1m", "1n7g", "1nax", "1oq5", "1osn", "1oyt", "1p2y", "1p5e", "1p62", "1t9b", "1tsl",
                        "1tt1", "1xdn", "2gpu", "2is7", "2iuz", "2p16", "2reg", "2w0s", "2yxj", "2zoz", "3byz", "3g2k",
                        "3i6O", "3O1H", "3pxf", "3r0t", "3shy", "3srr", "3t5y", "4alw", "4kya", "4pjt", "4qnb",
                        "4rao", "4rdl"]
    hajiebrahimi_ligand = ["GP6", "0QH", "E20", "D16", "115", "TQ3", "VIB", "AZM", "TNK", "STC", "BCZ",
                           "A3M", "NDP", "IH5", "CEL", "BVP", "FSN", "NCT", "TBS", "GEO", "1CS", "A15", "KAI", "ATP",
                           "OHT", "2CL", "D1H", "GG2", "CHT", "BVP", "N3C", "ET", "H11", "SKY", "GR6", "TMO", "2AN",
                           "FU9", "5FO", "Q20", "MCS", "HY7", "1UG", "2YQ", "1B0", "3L7", "FUC"]
    hajiebrahimi_chain = ["A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "H",
                          "A", "A", "B", "A", "A", "A", "A", "A", "A", "A", "A", "A", "B", "A", "B", "A", "A", "B", "B",
                          "A", "A", "A", "X", "B", "A", "E", "D", "A", "B", "A", ]

    hajiebrahimi = [list(combi) for combi in zip(hajiebrahimi_pdb, hajiebrahimi_ligand, hajiebrahimi_chain)]
    def run(element):
        print("Starting analysis of ", element[0])
        analysis = PdbIDAnalysis(element[0], element[1], data_path="../data/")
        return analysis.run(element[2], get_properties=False)

    haji_results = []
    for element in hajiebrahimi:
        result = run(element)
        if result is not None:
            haji_results.append(result)

    haji_results_columns = [list(entry.columns) for entry in haji_results]
    haji_proteins = [protein for protein in hajiebrahimi_pdb if (protein not in PdbIDAnalysis.reading_failed and
                                                                 protein not in PdbIDAnalysis.features_failed)]

    result_dataframe = pd.DataFrame.from_dict(dict(zip(haji_proteins, haji_results_columns)), orient="index")

    result_dataframe.to_pickle("results.pkl")

    analyze = StatisticalAnalysis("results.pkl")
    results_cleaned = analyze.get_results()
    baseline_cleaned = analyze.get_baseline(r"..\baseline.txt")
    pos, neg = analyze.compare(baseline=baseline_cleaned, results=results_cleaned)
