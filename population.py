'''this file is to design a Population class for evolutionary algorithm'''

from config import *
from molecule import Molecule

config = Config()

class Population(object):
    def __init__(self, ):
        self.population_pool = []
        self.population_size = config.poplution_size
        self.init_poplution_file_name = config.init_poplution_file_name
        self.crossover_rate = config.crossover_rate

        self._init_population()

    def _init_population(self):
        init_data = Chem.SDMolSupplier(self.init_poplution_file_name)
        for i in range(self.population_size):
            self.population_pool.append(Molecule(Chem.MolToSmiles(init_data[i])))

