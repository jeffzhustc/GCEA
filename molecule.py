'''This file is for unit of Graph based Chemical Evolutionary Algorithm(GCEA)'''
from config import Config
import rdkit.Chem as Chem
from rdkit.Chem.QED import qed
from utils import *
import numpy as np

#config = Config()

'''Molecule's object class be an individual in population for evolutionary algorithm'''


class Molecule(object):
    def __init__(self, smiles, config):
        self.smiles = smiles

        self.possible_bonds = config.possible_bonds
        self.table_of_elements = config.table_of_elements
        self.vocab_nodes_encode = config.vocab_nodes_encode
        self.mol = Chem.MolFromSmiles(smiles)

        self.adj = self._get_adj_mat(smiles)
        self.node_list = self._get_node_list(smiles)
        self.num_atom = len(self.node_list)
        self.expand_mat = self._get_expand_mat(self.adj, self.node_list)

        self.property = {
            'qed': qed(self.mol),
            'J_score': calc_score(self.mol)
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
            adj[[first], [second]] = self.possible_bonds.index(
                bond_list[first, second]) + 1
            adj[[second], [first]] = self.possible_bonds.index(
                bond_list[first, second]) + 1
        return adj

    def _get_node_list(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        G = mol2nx(mol)
        atomic_nums = nx.get_node_attributes(G, 'atomic_num')
        node_list = []
        for i in range(len(atomic_nums)):
            try:
                node_list.append(self.table_of_elements[atomic_nums[i]])
            except KeyError:
                pass
        return node_list

    def _get_expand_mat(self, adj, node_list):
        def _get_diag_mat(node_list):
            length = len(node_list)
            diag_mat = np.zeros([length, length])
            for i in range(length):
                diag_mat[[i], [i]] = self.vocab_nodes_encode[node_list[i]]
            return diag_mat

        diag_mat = _get_diag_mat(node_list)
        return adj + diag_mat
