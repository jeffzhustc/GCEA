'''this experiment is for comparison with GCEA's itself with different parameters'''

from molecule import *
from population import *
import time

import logging
import copy
from mutation import *
import rdkit.Chem as Chem
import rdkit.Chem.QED

import pickle

config = Config()
pop = Population(config)

parameters = {
    'pop_size': [100, 100, 500, 500],
    'is_mutate': [True, False, True, False],
    'logger_name': ['parameters'+str(i)+'.log' for i in range(4)],
    'time_count': 40
}


def get_logger(logger_name):
    logging.basicConfig(level=logging.INFO,
                        filename=logger_name,
                        datefmt='%Y/%m/%d %H:%M:%S',
                        format='%(asctime)s - %(name)s - %(levelname)s - %(lineno)d - %(module)s - %(message)s')
    logger = logging.getLogger(__name__)
    return logger


time_count = 0


def train(pop, logger, property_name='J_score', is_mutate=False, delta_time=120):
    time_count = 0
    starttime = (int)(time.time())

    for i in range(1000000):
        if time_count >= parameters['time_count']:
            break
        num_new_atom = pop.population_size
        temp_pool = []
        for i in range(num_new_atom):
            index_ = np.random.choice(num_new_atom, 2, replace=False)
            index1, index2 = index_[0], index_[1]
            mol1, mol2 = pop.population_pool[index1], pop.population_pool[index2]
            new_mol = pop.crossover(mol1, mol2)
            temp_pool.append(new_mol)

        if is_mutate:
            mutate_index = np.random.choice(num_new_atom, int(
                num_new_atom*config.mutate_rate), replace=False)
            for i in mutate_index:
                temp_molecule = copy.deepcopy(pop.population_pool[i])
                new_molecule = mutate(temp_molecule)
                if temp_molecule != new_molecule:
                    temp_pool.append(new_molecule)

        total = temp_pool + pop.population_pool
        total = set(total)
        res = sorted(
            total, key=lambda x: x.property[property_name], reverse=True)

        res_score = [res[i].property[property_name]
                     for i in range(pop.population_size)]
        print('mean is : ' + str(np.mean(res_score)) + ' std is : ' +
              str(np.std(res_score)) + ' max is : ' + str(np.max(res_score)))

        pop.population_pool = []
        for i in range(pop.population_size):
            pop.population_pool.append(res[i])
        # max_pool.append(res[0])

        currenttime = (int)(time.time())
        if (currenttime - starttime) > delta_time:
            starttime = currenttime
            print('ok')
            logger.info('mean is : ' + str(np.mean(res_score)) + ' std is : ' +
                        str(np.std(res_score)) + ' max is : ' + str(np.max(res_score)))
            time_count += 1
            if np.std(res_score) < 1e-6:
                break


for i in range(4):
    config.poplution_size = parameters['pop_size'][i]
    pop = Population(config)
    logger = get_logger(parameters['logger_name'][i])
    is_mutate = parameters['is_mutate'][i]

    train(pop, logger, is_mutate=is_mutate)
