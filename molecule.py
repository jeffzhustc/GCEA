'''This file is for unit of Graph based Chemical Evolutionary Algorithm(GCEA)'''
from config import Config
import rdkit.Chem as Chem
from rdkit.Chem.QED import qed
from score_util import calc_score
from utils import *
import numpy as np

config = Config()

class Molecule(object):
    def __init__(self, smiles):
        self.smiles = smiles

        self.adj = self._get_adj_mat(smiles)
        self.node_list = self._get_node_list(smiles)
        self.num_atom = len(self.node_list)
        self.diag_mat = self._get_diag_mat(self.node_list)
        self.expand_mat = self._get_expand_mat(self.adj, self.diag_mat)
        self.mol = Chem.MolFromSmiles(smiles)

        self.property = {
            'qed' : qed(self.mol),
            'J_socre' : calc_score(self.mol)
        }

    def _get_adj_mat(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        Chem.Kekulize(mol)
        G = mol2nx(mol)
        atomic_nums = nx.get_node_attributes(G, 'atomic_num')
        adj = np.zeros([len(atomic_nums), len(atomic_nums)])
        bond_list = nx.get_edge_attributes(G, 'bond_type')
        for edge in G.edges():
            first, second = edge
            adj[[first], [second]] = possible_bonds.index(bond_list[first, second]) + 1
            adj[[second], [first]] = possible_bonds.index(bond_list[first, second]) + 1
        return adj