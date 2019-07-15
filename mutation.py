'''this file is for mutation operation in evolutionary algorithm'''

from utils import *
from molecule import Molecule

deconfig = Config()

def _add_bond(molecule):
    if molecule.num_atom < 2:
        return molecule

    temp_expand_adj = molecule.expand_mat
    temp_adj = molecule.adj
    mask_row = mask(temp_expand_adj)

    goal_mol = None

    for i in mask_row:
        row = temp_adj[i]
        for j in range(len(row)):
            if row[j] > 0 and j in mask_row:
                temp_adj[i][j] += 1
                temp_adj[j][i] += 1
                goal_adj = temp_adj
                goal_node_list = molecule.node_list
                goal_mol = adj2mol(goal_node_list, goal_adj.astype(int), deconfig.possible_bonds)
                goal_smiles = Chem.MolToSmiles(goal_mol)
                break
        if goal_mol != None:
            break

    if goal_mol != None:
        return Molecule(goal_smiles)
    else:
        return molecule