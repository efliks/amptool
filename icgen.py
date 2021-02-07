from collections import namedtuple

from rdkit import Chem
from rdkit.Chem import rdMolTransforms


IntCoor = namedtuple('IntCoor', [
    'i', 'j', 'k', 'l', 'ij', 'ijk', 'ijkl'
])


def sort_heavy_atoms_r(mol, idx_list):
    last_atom_idx = idx_list[-1]
    last_atom = mol.GetAtomWithIdx(last_atom_idx)
    for neighbor in last_atom.GetNeighbors():
        if neighbor.GetSymbol() != 'H':
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in idx_list:
                idx_list.append(neighbor_idx)
                sort_heavy_atoms_r(mol, idx_list)


def add_hydrogens(mol, idx_list):
    for atom_idx in idx_list:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() == 'H':
                neighbor_idx = neighbor.GetIdx()
                idx_list.append(neighbor_idx)


def get_predecessors_r(mol, atom_idx, order_map, predecessor_list):
    atom = mol.GetAtomWithIdx(atom_idx)
    for neighbor in atom.GetNeighbors():
        neighbor_idx = neighbor.GetIdx()
        if order_map[neighbor_idx] < order_map[atom_idx]:
            predecessor_list.append(neighbor_idx)
            if len(predecessor_list) == 3:
                return True
            is_done = get_predecessors_r(
                mol, neighbor_idx, order_map, predecessor_list
            )
            if is_done:
                return True
    return False


def make_internal_map(mol, idx_list, order_map):
    internal_map = {}
    for atom_idx in idx_list:
        preceding_atoms = []
        get_predecessors_r(mol, atom_idx, order_map, preceding_atoms)
        internal_map[atom_idx] = preceding_atoms
    return internal_map


def calculate_internal(mol, idx_list, internal_map):
    internal_coords = []
    conformer = mol.GetConformer(0)
    for atom_idx in idx_list:
        other_atom_idxs = internal_map[atom_idx]
        if len(other_atom_idxs) == 3:
            ij = atom_idx, other_atom_idxs[0]
            ijk = atom_idx, other_atom_idxs[0], other_atom_idxs[1]
            ijkl = atom_idx, other_atom_idxs[0], other_atom_idxs[1], other_atom_idxs[2]
            bond = rdMolTransforms.GetBondLength(conformer, *ij)
            angle = rdMolTransforms.GetAngleDeg(conformer, *ijk)
            torsion = rdMolTransforms.GetDihedralDeg(conformer, *ijkl)
            internal_coords.append(IntCoor(
                atom_idx, *other_atom_idxs, bond, angle, torsion
            ))
    return internal_coords


mol_h = Chem.MolFromPDBFile('monomer.pdb', removeHs=False)

label_idx_map = {}
for atom in mol_h.GetAtoms():
    pdb_info = atom.GetPDBResidueInfo()
    atom_label = pdb_info.GetName().strip()
    atom_idx = atom.GetIdx()
    label_idx_map[atom_label] = atom_idx

order_list = [label_idx_map['O']]
sort_heavy_atoms_r(mol_h, order_list)
add_hydrogens(mol_h, order_list)

idx_label_map = {idx: label for label, idx in label_idx_map.items()}
order_list_labels = list(map(idx_label_map.get, order_list))

mol_h_reorder = Chem.RenumberAtoms(mol_h, order_list)
Chem.MolToPDBFile(mol_h_reorder, 'monomer_reorder.pdb')

idx_order_map = {
    atom_idx: order_num
    for order_num, atom_idx in enumerate(order_list, 1)
}
internal_map = make_internal_map(mol_h, order_list, idx_order_map)
internal_coords = calculate_internal(mol_h, order_list, internal_map)

labeled_internal_coords = []
for internal_coord in internal_coords:
    atom_labels = tuple(map(idx_label_map.get, (
        internal_coord.i, internal_coord.j, internal_coord.k, internal_coord.l
    )))
    labeled_internal_coords.append(IntCoor(
        *atom_labels, internal_coord.ij, internal_coord.ijk, internal_coord.ijkl
    ))

for crd in labeled_internal_coords:
    print(f'{crd.i} {crd.j} {crd.k} {crd.l} {crd.ij} {crd.ijk} {crd.ijkl}')
