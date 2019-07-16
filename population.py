'''this file is to design a Population class for evolutionary algorithm'''
from config import *
from molecule import Molecule
import numpy as np
from queue import Queue
from utils import *
from random import *

class Population(object):
    def __init__(self, config):
        self.population_pool = []
        self.population_size = config.poplution_size
        self.init_poplution_file_name = config.init_poplution_file_name
        self.crossover_rate = config.crossover_rate
        self.vocab_nodes_encode = config.vocab_nodes_encode
        self.crossover_mu = config.crossover_mu
        self.crossover_sigma = config.crossover_sigma
        self.graph_size = config.graph_size
        self.possible_bonds = config.possible_bonds

        self._init_population()

    # init population using 500 random molecules from zinc dataset
    def _init_population(self):
        init_data = Chem.SDMolSupplier(self.init_poplution_file_name)
        for i in range(self.population_size):
            self.population_pool.append(Molecule(Chem.MolToSmiles(init_data[i]), config))

    # using BFS to sample a subgraph from parent molecule
    def _bfs_molecule(self, molecule):
        length = molecule.num_atom

        sigma = self.crossover_sigma * length
        mu = self.crossover_mu * length
        '''using normal distribution to select the number of molecules'''
        '''this mu and sigma of normal distribution is decided by the above paramters'''
        num_sample = (int)(np.random.normal(loc=mu, scale=sigma))
        while num_sample >= length:
            num_sample = (int)(np.random.normal(loc=mu, scale=sigma))

        '''we choose the start node from the set who has the most degrees for BFS algorithm '''
        adj = molecule.adj
        np.putmask(adj, adj >= 1, 1)
        row_sum = np.sum(adj, axis=0)
        max_ = np.max(row_sum)
        index_arr = [i for i in range(len(row_sum)) if row_sum[i] == max_]
        index_start = np.random.choice(index_arr, 1)[0]

        hash = np.zeros(length)
        res = []
        q = Queue(self.graph_size)
        q.put(index_start)
        res.append(index_start)
        hash[index_start] = 1

        while len(res) < num_sample:
            node = q.get()
            node_list = list(np.squeeze(np.argwhere(adj[node] >= 1), axis=1))
            for n in node_list:
                if hash[n] != 1:
                    hash[n] = 1
                    q.put(n)
                    res.append(n)

        temp_node_list = []
        '''res store the index of atom in the molecules'''
        for i in res:
            temp_node_list.append(molecule.node_list[i])

        '''reconstruct the subgraph in adj form'''
        num_atom = len(res)
        temp_mat = np.zeros([num_atom, num_atom])
        #print(res)
        for i in range(num_atom - 1):
            for j in range(1, num_atom):
                temp_mat[i, j] = molecule.expand_mat[res[i], res[j]]
                temp_mat[j, i] = molecule.expand_mat[res[j], res[i]]

        temp_mat = temp_mat.astype(int)
        length = len(temp_node_list)

        for i in range(length):
            temp_mat[i,i] = self.vocab_nodes_encode[temp_node_list[i]]
        # mol = adj2mol(temp_node_list, temp_mat, possible_bonds)
        # Draw.MolToFile(mol, 'test2.png')
        '''temp_mat for expand_mat, temp_node_list for node_list'''
        return temp_mat, temp_node_list

    # combine two subgraph to get an new graph represented as a molecule
    def _combine_two_subgraph(self, subg1, subg2):
        len1, len2 = subg1.shape[0], subg2.shape[0]
        goal_adj = np.zeros([len1 + len2, len1 + len2])
        goal_adj[:len1, :len1] = subg1
        goal_adj[len1:, len1:] = subg2
        ''' this part can be modified '''
        row1, row2 = mask(subg1), mask(subg2)
        row2 = [i + len1 for i in row2]
        index1, index2 = choice(row1), choice(row2)
        goal_adj[index1, index2] = 1
        goal_adj[index2, index1] = 1
        return goal_adj

    # we sample from two molecule and get an new molecule after mutation
    def crossover(self, mol1, mol2):
        temp_mat1, temp_node_list1 = self._bfs_molecule(mol1)
        temp_mat2, temp_node_list2 = self._bfs_molecule(mol2)

        goal_mat = self._combine_two_subgraph(temp_mat1, temp_mat2).astype(int)
        goal_list = temp_node_list1 + temp_node_list2
        adj = goal_mat
        for i in range(len(goal_list)):
            adj[i,i] = 0
        mol_temp = adj2mol(goal_list, adj, self.possible_bonds)

        return Molecule(Chem.MolToSmiles(mol_temp), config)


# import rdkit.Chem.Draw as Draw
#
# config = Config()
# pop_test = Population(config)
# mol_test = pop_test.population_pool[0]
# test_mat, test_node_list = pop_test._bfs_molecule(mol_test)
# new_mol = adj2mol(test_node_list, test_mat, config.possible_bonds)
# Draw.MolToFile(mol_test.mol, 'old_mol.png')
# Draw.MolToFile(new_mol, 'new_mol.png')