'''this experiment is for the comparision with ChemGE'''

from molecule import *
from population import *
import time

import logging
import rdkit.Chem as Chem
import rdkit.Chem.QED

import pickle

logging.basicConfig(level=logging.INFO,
                    filename='experiment_add_smi_vae.log',
                    datefmt='%Y/%m/%d %H:%M:%S',
                    format='%(asctime)s - %(name)s - %(levelname)s - %(lineno)d - %(module)s - %(message)s')
logger = logging.getLogger(__name__)

config = Config()
pop = Population(config)

max_pool = []

def train(pop, property_name='J_score', is_mutate=False):
    # res_qed = [pop.population_pool[i].J_score for i in range(pop.population_size)]

    starttime = (int)(time.time())

    for i in range(100):
        num_new_atom = pop.population_size
        temp_pool = []
        for i in range(num_new_atom):
            index_ = np.random.choice(num_new_atom, 2, replace=False)
            index1, index2 = index_[0], index_[1]
            #print(str(index1) + ' ' + str(index2))
            mol1, mol2 = pop.population_pool[index1], pop.population_pool[index2]
            new_mol = pop.crossover(mol1, mol2)
            temp_pool.append(new_mol)

        if is_mutate:
            mutate_index = np.random.choice(num_new_atom, int(num_new_atom*config.mutate_rate), replace=False)

        total = temp_pool + pop.population_pool
        total = set(total)
        res = sorted(total, key=lambda x:x.property[property_name], reverse=True)

        res_score = [res[i].property[property_name] for i in range(pop.population_size)]
        print('mean is : '+ str(np.mean(res_score)) + ' std is : ' + str(np.std(res_score)) + ' max is : '+ str(np.max(res_score)))

        pop.population_pool = []
        for i in range(pop.population_size):
            pop.population_pool.append(res[i])
        max_pool.append(res[0])

        currenttime = (int)(time.time())
        if (currenttime - starttime) > 1800:
            starttime = currenttime
            logger.info('mean is : '+ str(np.mean(res_score)) + ' std is : '+ str(np.std(res_score)) + ' max is : ' + str(np.max(res_score)))

train(pop)
