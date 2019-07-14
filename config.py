'''This file is for hypeparameters and configuration'''
import rdkit.Chem as Chem

class Config(object):
    def __init__(self):
        self.possible_bonds = [
            Chem.rdchem.BondType.SINGLE,
            Chem.rdchem.BondType.DOUBLE,
            Chem.rdchem.BondType.TRIPLE
        ]