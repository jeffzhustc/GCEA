from score_util import *
import numpy as np
from utils import *
from molecule import Molecule

deconfig = Config()

def mutate(molecule):
    mutate_rate = deconfig.mutation_rate
    choice = np.random.choice(3, 1, mutate_rate)
    if choice[0] == 0:
        return _add_atom(molecule)
    if choice[0] == 1:
        return _add_atom_between_bond(molecule)
    if choice[0] == 2:
        return _add_bond(molecule)

def _add_bond(molecule):
    if molecule.num_atom < 2:
        return molecule

    temp_expand_adj = molecule.expand_mat
    temp_adj = molecule.adj
    mask_row = mask(temp_expand_adj)

    goal_mol = None
    goal_smiles = None

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
        return Molecule(goal_smiles, deconfig)
    else:
        return molecule

def _add_atom_between_bond(molecule):
    temp_elements = deconfig.temp_elements
    atom_index = np.random.choice(deconfig.length_elements, 1)[0]
    atom = temp_elements[atom_index]

    temp_adj = molecule.adj

    length = molecule.num_atom
    insert_index1 = np.random.choice(length, 1)
    insert_row = temp_adj[insert_index1][0]

    insert_index2 = 0
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

    return Molecule(goal_smiles, deconfig)

def _add_atom(molecule):
    if molecule.num_atom < 1:
        return molecule

    temp_node_list = molecule.node_list
    #print(len(temp_node_list))
    temp_expand_adj = molecule.expand_mat
    #print(temp_expand_adj.shape[0])
    temp_adj = molecule.adj

    temp_elements = deconfig.temp_elements

    atom_index = np.random.choice(deconfig.length_elements, 1)[0]
    atom = temp_elements[atom_index]
    mask_row = mask(temp_expand_adj)
    if len(mask_row) < 1:
        return molecule
    mask_index = np.random.choice(mask_row, 1)[0]

    goal_length = molecule.num_atom + 1
    goal_adj = np.zeros([goal_length, goal_length])
    goal_adj[:goal_length-1 , :goal_length-1] = temp_adj
    goal_adj[goal_length-1, mask_index] = goal_adj[mask_index, goal_length-1] = 1

    temp_node_list.append(atom)
    goal_node_list = temp_node_list

    goal_mol = adj2mol(goal_node_list, goal_adj.astype(int), deconfig.possible_bonds)
    goal_smiles = Chem.MolToSmiles(goal_mol)
    return Molecule(goal_smiles, deconfig)