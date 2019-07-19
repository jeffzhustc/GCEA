import pickle
import random
import rdkit.Chem as Chem

def change_random_molecule(config):
    random_size = config.poplution_size
    data = pickle.load(open('/home/jeffzhu/aaai_ga/data/zinc_clean_smi.pkl', 'rb'))
    length = len(data)

    random_index = random.choice(length, random_size)

    random_pool = [data[i] for i in random_index]
    w = Chem.SDWriter('/home/jeffzhu/aaai_ga/data/random.sdf')

    for i in random_pool:
        w.write(i)

        