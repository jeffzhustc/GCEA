'''This file is for hypeparameters and configuration'''
import rdkit.Chem as Chem

class Config(object):
    def __init__(self):
        '''hyperparameters for Molecule class'''
        self.possible_bonds = [
            Chem.rdchem.BondType.SINGLE,
            Chem.rdchem.BondType.DOUBLE,
            Chem.rdchem.BondType.TRIPLE
        ]
        self.table_of_elements = {
            6 : 'C',
            7 : 'N',
            8 : 'O',
            9 : 'F',
            16 : 'S',
            17 : 'Cl',
            35 : 'Br',
            53 : 'I',
        }
        self.vocab_nodes_encode = {
            'C': 1, 'N': 2, 'O': 3, 'S': 3, 'F': 4, 'Cl': 4, 'Br': 4, 'I': 4
        }
        '''hyperparameters for Population class'''
        self.poplution_size = 500
        self.crossover_rate = 0.8
        self.init_poplution_file_name = 'data/init_population.sdf'