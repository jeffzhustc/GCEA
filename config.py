'''This file is for hypeparameters and configuration'''
import rdkit.Chem as Chem

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