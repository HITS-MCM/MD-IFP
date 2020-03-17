import re
import wget
import os
from Bio.PDB import *


class AutomatedAnalysis:
    def __init__(self, data_path='../Data/', verbose=False):
        self.data_path = data_path
        self.verbose = verbose

# Method downloads the given pdb dataset from the RCSB. It checks for correct formatting and folder structure
    def download_pdb(self, wanted_data):
        if (re.match(r"[1-9]\w\w\w", wanted_data) is None) or len(wanted_data) != 4:
            raise ValueError('Input data does not conform to pdb standard')
        save_path = self.data_path + wanted_data + '.pdb'
        if not os.path.isdir(self.data_path):
            os.mkdir(self.data_path)
        if os.path.isfile(save_path) is True:
            print('File {} has already been downloaded.'.format(wanted_data))
        else:
            if self.verbose is True:
                print('Downloading {}.pdb from RCSB to {}'.format(wanted_data, save_path))
            wget.download('https://files.rcsb.org/download/{0}.pdb'.format(wanted_data), out=save_path)

# TODO: Add commentary and explanation
    def split_pdb_automated(self, pdb_file, pdb_chain, pdb_ligand):
        structure = PDBParser().get_structure(pdb_file, self.data_path + pdb_file + '.pdb')

        class LigandSelect(Select):
            def accept_residue(self, residue):
                if residue.get_resname() == pdb_ligand:
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
        split_path = self.data_path + pdb_file + "split/"
        io = PDBIO()
        io.set_structure(structure)
        if self.verbose is True:
            print("Extracting protein from residue {0} to {1}, ligand and water. Saving to {2}".format(
                1, len(structure[0][pdb_chain]), split_path))
        if not os.path.isdir(split_path):
            os.mkdir(split_path)
        io.save(split_path + "protein" + '.pdb', Dice.ChainSelector("A", 1, len(structure[0][pdb_chain])))
        io.save(split_path + pdb_ligand + '.pdb', LigandSelect())
        io.save(split_path + 'water' + '.pdb', WaterSelect())

# TODO: Should we keep this or can we remove it?
    def split_pdb_manual(self, pdb_file):
        with open(self.data_path + pdb_file + '.pdb') as p_file:
            protein = ''
            ligand = ''
            for line in p_file:
                re1_protein = bool(re.match(r'(^\bATOM)', line))
                re2_protein = bool(re.match(r'ATOM[1-9\s]', line))
                re1_ligand = bool(re.match(r'(^\bHETATM)', line))
                re2_ligand = bool(re.match(r'HETATM[1-9\s]', line))
                if re1_protein and re2_protein:
                    protein = protein + line
                    continue
                elif re1_ligand and re2_ligand:
                    ligand = ligand + line
        split_path = self.data_path+pdb_file+'split'
        if not os.path.isdir(split_path):
            os.mkdir(split_path)
        with open(split_path + '/protein.pdb', mode='a+') as protein_pdb:
            protein_pdb.write(protein)
        with open(split_path + '/ligand.pdb', mode='a+') as ligand_pdb:
            ligand_pdb.write(ligand)


if __name__ == '__main__':
    analyze = AutomatedAnalysis(verbose=True)
    analyze.download_pdb("2hcp")
    analyze.download_pdb("3hbp")
    analyze.split_pdb_automated("2hcp", "A", "BIC")
    analyze.split_pdb_automated("3hbp", "A", "BIC")
