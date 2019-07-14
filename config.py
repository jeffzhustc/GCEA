'''This file is for hypeparameters and configuration'''
import rdkit.Chem as Chem

'''restore paramters using a class'''
class Config(object):
    def __init__(self):
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
        self.poplution_size = 500
        self.crossover_rate = 0.8
        self.init_poplution_file_name = 'data/init_population.sdf'
        self.crossover_mu = 0.5
        self.crossover_sigma = 0.1
        self.graph_size = 80
        self.full_valence = 5