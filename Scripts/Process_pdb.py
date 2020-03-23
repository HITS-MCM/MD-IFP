import re
import wget
import os
import subprocess
import MDAnalysis as mda
from Bio.PDB import *
from Scripts.IFP_generation import IFP, table_combine
from rdkit import Chem, RDConfig
from rdkit.Chem.rdMolChemicalFeatures import MolChemicalFeatureFactory


class AutomatedAnalysis:
    def __init__(self, protein, ligand, data_path='../Data/', verbose=False):
        self.data_path = data_path
        self.verbose = verbose
        self.protein = protein
        self.ligand = ligand

# Method downloads the given pdb dataset from the RCSB. It checks for correct formatting and folder structure
    def download_pdb(self):
        if (re.match(r"[1-9]\w\w\w", self.protein) is None) or len(self.protein) != 4:
            raise ValueError('Input data does not conform to pdb standard')
        save_path = self.data_path + self.protein + '.pdb'
        if not os.path.isdir(self.data_path):
            os.mkdir(self.data_path)
        if os.path.isfile(save_path) is True:
            print('File {} has already been downloaded.'.format(self.protein))
        else:
            if self.verbose is True:
                print('Downloading {}.pdb from RCSB to {}'.format(self.protein, save_path))
            wget.download('https://files.rcsb.org/download/{0}.pdb'.format(self.protein), out=save_path)

# TODO: Add commentary and explanation
    def split_pdb_automated(self, pdb_chain):
        structure = PDBParser().get_structure(self.protein, self.data_path + self.protein + '.pdb')
        ligand = self.ligand

        class LigandSelect(Select):
            def accept_residue(self, residue):
                if residue.get_resname() == ligand:
                    return 1
                else:
                    return 0

            def accept_chain(self, chain):
                if chain.get_id() == pdb_chain:
                    return 1
                else:
                    return 0

        class WaterSelect(Select):
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
        split_path = self.data_path + self.protein + "split/"
        io = PDBIO()
        io.set_structure(structure)
        if self.verbose is True:
            print("Extracting protein from residue {0} to {1}, ligand and water. Saving to {2}".format(
                1, len(structure[0][pdb_chain]), split_path))
        if not os.path.isdir(split_path):
            os.mkdir(split_path)
        io.save(split_path + "protein" + '.pdb', Dice.ChainSelector("A", 1, len(structure[0][pdb_chain])))
        io.save(split_path + self.ligand + '.pdb', LigandSelect())
        io.save(split_path + 'water' + '.pdb', WaterSelect())

    def prepare_data(self, chimera_path="C:/Program Files/Chimera 1.14/bin/chimera.exe"):
        subprocess.run("\"{0}\" --nogui --script \"{1}/chimera_hydrogen_mol2.py {2}/{3}split/\"".format(
            chimera_path, os.getcwd(), self.data_path, self.protein), shell=True)
        with open(self.data_path + self.protein + "split/merged.pdb", "w+") as outfile:
            for file in ["protein", self.ligand, "water"]:
                with open(self.data_path + self.protein + "split/" + file + ".pdb", "r") as infile:
                    lines = [line for line in infile if line.split()[0] not in ["CONECT", "MASTER", "END", "TER"]]
                    outfile.writelines(lines)

    def calculate_ifp(self, get_properties=False):
        ligand_pdb = self.data_path + self.protein + "split/" + self.ligand + ".pdb"
        ligand_mol2 = self.data_path + self.protein + "split/" + self.ligand + ".mol2"
        pdb_merged = self.data_path + self.protein + "split/merged.pdb"
        with open(ligand_pdb, 'r') as fasta_file:
            list_labels = [line.split()[2] for line in fasta_file.readlines() if line.split()[0] in ["ATOM", "HETATM"]]
        mol = Chem.rdmolfiles.MolFromMol2File(ligand_mol2, removeHs=False)
        factory_name = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
        factory = Chem.ChemicalFeatures.BuildFeatureFactory(factory_name)
        features = factory.GetFeaturesForMol(mol)
        properties = {feature.GetFamily(): [] for feature in features}
        for feature in features:
            current_atoms = [list_labels[atom_id] for atom_id in list(feature.GetAtomIds())]
            for atom in current_atoms:
                properties[feature.GetFamily()].append(atom)
        u = mda.Universe(pdb_merged)
        df_prop, df_hb, df_wb = IFP(u, self.ligand, properties, RE=False)
        df_prop_complete = table_combine(df_hb, df_wb, df_prop, self.ligand)
        if get_properties is True:
            return df_prop_complete, properties
        return df_prop_complete


if __name__ == '__main__':
    analyze = AutomatedAnalysis("6ei5", "B5Q", verbose=True)
    analyze.download_pdb()
    analyze.split_pdb_automated("A")
    analyze.prepare_data()
    test = analyze.calculate_ifp()
