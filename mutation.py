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

def _add_atom_between_bond(molecule):
    temp_elements = {
        0: 'C',
        1: 'N',
        2: 'O',
        3: 'S'
    }
    atom_index = np.random.choice(4, 1)[0]
    atom = temp_elements[atom_index]

    temp_adj = molecule.adj

    length = molecule.num_atom
    insert_index1 = np.random.choice(length, 1)
    insert_row = temp_adj[insert_index1][0]

    for i in range(len(insert_row)):
        if insert_row[i] > 0:
            insert_index2 = i

    temp_adj[insert_index1, insert_index2] = temp_adj[insert_index2, insert_index1] = 0

    goal_adj = np.zeros([length+1, length+1])
    goal_adj[:length, :length] = temp_adj
    goal_adj[length, insert_index1] = goal_adj[insert_index1, length] = 1
    goal_adj[insert_index2, length] = goal_adj[length, insert_index2] = 1

    temp_node_list = molecule.node_list
    temp_node_list.append(atom)
    goal_node_list = temp_node_list

    goal_mol = adj2mol(goal_node_list, goal_adj.astype(int), deconfig.possible_bonds)
    goal_smiles = Chem.MolToSmiles(goal_mol)

    return Molecule(goal_smiles)