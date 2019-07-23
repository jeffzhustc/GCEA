from rdkit import Chem
from rdkit.Chem import Draw


def Draw_smi(smi, file_name):
    mol = Chem.MolFromSmiles(smi)
    img = Draw.MolToFile(mol, file_name)

    return img
