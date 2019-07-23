'''this file is just for general function'''
import rdkit.Chem as Chem
import networkx as nx
from config import Config
from rdkit.Chem import Descriptors
from rdkit.Chem import MolFromSmiles
from rdkit.Chem import rdmolops
from rdkit import rdBase

import sascorer
import numpy as np

config = Config()

# mask the matrix to find valid row to take mutation or crossover operation


def mask(expand_adj):
    node_num = np.count_nonzero(expand_adj.diagonal())
    row_sum = np.sum(expand_adj[:node_num, :node_num], axis=0)
    mask_row = np.argwhere(
        row_sum < config.full_valence).squeeze(axis=1).tolist()
    return mask_row

# adj2mol is to convert adjacent matrix into mol object in rdkit


def adj2mol(nodes, adj, possible_bonds):
    mol = Chem.RWMol()

    for i in range(len(nodes)):
        # print(nodes[i])
        atom = Chem.Atom(nodes[i])
        mol.AddAtom(atom)

    for i in range(len(nodes)-1):
        for j in range(i + 1, len(nodes)):
            if adj[i, j]:
                mol.AddBond(i, j, possible_bonds[adj[i, j] - 1])

    return mol

# mol2nx is to convert mol object in rdkit into network object


def mol2nx(mol):
    G = nx.Graph()

    for atom in mol.GetAtoms():
        G.add_node(atom.GetIdx(),
                   atomic_num=atom.GetAtomicNum(),
                   formal_charge=atom.GetFormalCharge(),
                   chiral_tag=atom.GetChiralTag(),
                   hybridization=atom.GetHybridization(),
                   num_explicit_hs=atom.GetNumExplicitHs(),
                   is_aromatic=atom.GetIsAromatic())
    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(),
                   bond.GetEndAtomIdx(),
                   bond_type=bond.GetBondType())
    return G


def calc_score(mol):
    logP_mean = 2.457  # np.mean(logP_values)
    logP_std = 1.434  # np.std(logP_values)
    SA_mean = -3.053  # np.mean(SA_scores)
    SA_std = 0.834  # np.std(SA_scores)
    cycle_mean = -0.048  # np.mean(cycle_scores)
    cycle_std = 0.287  # np.std(cycle_scores)

    molecule = mol
    if Descriptors.MolWt(molecule) > 500:
        return -1e10
    current_log_P_value = Descriptors.MolLogP(molecule)
    current_SA_score = -sascorer.calculateScore(molecule)
    cycle_list = nx.cycle_basis(
        nx.Graph(rdmolops.GetAdjacencyMatrix(molecule)))
    if len(cycle_list) == 0:
        cycle_length = 0
    else:
        cycle_length = max([len(j) for j in cycle_list])
    if cycle_length <= 6:
        cycle_length = 0
    else:
        cycle_length = cycle_length - 6
    current_cycle_score = -cycle_length

    current_SA_score_normalized = (current_SA_score - SA_mean) / SA_std
    current_log_P_value_normalized = (
        current_log_P_value - logP_mean) / logP_std
    current_cycle_score_normalized = (
        current_cycle_score - cycle_mean) / cycle_std

    score = (current_SA_score_normalized
             + current_log_P_value_normalized
             + current_cycle_score_normalized)