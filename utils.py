'''this file is just for general function'''
import rdkit.Chem as Chem
import networkx as nx
from config import Config
import numpy as np

config = Config()

# mask the matrix to find valid row to take mutation or crossover operation
def mask(expand_adj):
    node_num = np.count_nonzero(expand_adj.diagonal())
    row_sum = np.sum(expand_adj[:node_num, :node_num], axis=0)
    mask_row = np.argwhere(row_sum < config.full_valence).squeeze(axis=1).tolist()
    return mask_row

# adj2mol is to convert adjacent matrix into mol object in rdkit
def adj2mol(nodes, adj, possible_bonds):
    mol = Chem.RWMol()

    for i in range(len(nodes)):
        #print(nodes[i])
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