'''this file is just for '''
import rdkit.Chem as Chem

# adj2mol is to convert adjacent matrix into mol object in rdkit
def adj2mol(nodes, adj, possible_bonds):
    mol = Chem.RWMol()

    for i in range(len(nodes)):
        #print(nodes[i])
        atom = Chem.Atom(nodes[i])
        mol.AddAtom(atom)

    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            if adj[i, j]:
                mol.AddBond(i, j, possible_bonds[adj[i, j] - 1])

    return mol

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