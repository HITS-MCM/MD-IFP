import re
import wget
import os
import Bio


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

    def split_pdb(self, pdb_file):
        with open(self.data_path + pdb_file + '.pdb') as p_file:
            protein = ''
            ligand = ''
            for line in p_file:
                re1_protein = bool(re.match(r'(^\bATOM)', line))
                re2_protein = bool(re.match(r'ATOM[1-9\s]', line))
                re1_ligand = bool(re.match(r'(^\bHETATOM)', line))
                re2_ligand = bool(re.match(r'HETATOM[1-9\s]', line))
                if re1_protein and re2_protein:
                    protein = protein + line
                    continue
                elif re1_ligand and re2_ligand:
                    ligand = ligand + line
        split_path = self.data_path+pdb_file+'split'
        if not os.path.isdir(split_path):
            os.mkdir(split_path)
        with open(split_path + 'protein.txt', mode='a+') as protein_pdb:
            protein_pdb.write(protein)
        with open(split_path + 'ligand.txt', mode='a+') as ligand_pdb:
            ligand_pdb.write(ligand)






if __name__ == '__main__':
    analyze = AutomatedAnalysis(verbose=True)
    analyze.download_pdb('2hcp')
    analyze.download_pdb('3hpb')
    analyze.split_pdb('test')
