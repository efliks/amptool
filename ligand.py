import json
import networkx as nx
from collections import namedtuple, defaultdict

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import BondType


CharmmAtom = namedtuple('CharmmAtom', [
    'label',
    'ff_type',
    'charge',
    'charge_group_id',
    'serial'
])

atom_map = {
   17: CharmmAtom('N'  , 'NH1',  -0.47 ,   1,  -1),
   48: CharmmAtom('H'  , 'H'  ,   0.15 ,   1,  -1),
   49: CharmmAtom('H1' , 'H'  ,   0.16 ,   1,  -1),   # added proton
   16: CharmmAtom('CA' , 'CT1',   0.07 ,   1,  -1),
   47: CharmmAtom('HA' , 'HB' ,   0.09 ,   1,  -1),   # end group
   15: CharmmAtom('CB' , 'CT2',  -0.18 ,   2,  -1),
   45: CharmmAtom('HB1', 'HA' ,   0.09 ,   2,  -1),
   46: CharmmAtom('HB2', 'HA' ,   0.09 ,   2,  -1),   # end group
   14: CharmmAtom('CG' , 'CT2',  -0.28 ,   3,  -1),
   43: CharmmAtom('HG1', 'HA' ,   0.09 ,   3,  -1),
   44: CharmmAtom('HG2', 'HA' ,   0.09 ,   3,  -1),
   12: CharmmAtom('CD' , 'CC' ,   0.62 ,   3,  -1),
   13: CharmmAtom('OE1', 'OC' ,  -0.76 ,   3,  -1),
   11: CharmmAtom('OE2', 'OC' ,  -0.76 ,   3,  -1),   # end group
   18: CharmmAtom('C'  , 'C'  ,   0.51 ,   4,  -1),
   20: CharmmAtom('O'  , 'OB' ,  -0.51 ,   4,  -1),   # added atom
   19: CharmmAtom('O1' , 'OH1',  -0.61 ,   4,  -1),   # added atom
   50: CharmmAtom('H11', 'H'  ,   0.61 ,   4,  -1),   # end group
   10: CharmmAtom('CI' , 'CT2',  -0.18 ,   5,  -1),
   41: CharmmAtom('HI1', 'HA' ,   0.09 ,   5,  -1),
   42: CharmmAtom('HI2', 'HA' ,   0.09 ,   5,  -1),   # end group
    9: CharmmAtom('CJ' , 'CT2',  -0.18 ,   6,  -1),
   39: CharmmAtom('HJ1', 'HA' ,   0.09 ,   6,  -1),
   40: CharmmAtom('HJ2', 'HA' ,   0.09 ,   6,  -1),   # end group
    8: CharmmAtom('CK' , 'CT2',  -0.18 ,   7,  -1),
   37: CharmmAtom('HK1', 'HA' ,   0.09 ,   7,  -1),
   38: CharmmAtom('HK2', 'HA' ,   0.09 ,   7,  -1),   # end group
    7: CharmmAtom('CL' , 'CT2',  -0.18 ,   8,  -1),
   35: CharmmAtom('HL1', 'HA' ,   0.09 ,   8,  -1),
   36: CharmmAtom('HL2', 'HA' ,   0.09 ,   8,  -1),   # end group
    6: CharmmAtom('CM' , 'CT2',  -0.18 ,   9,  -1),
   33: CharmmAtom('HM1', 'HA' ,   0.09 ,   9,  -1),
   34: CharmmAtom('HM2', 'HA' ,   0.09 ,   9,  -1),   # end group
    5: CharmmAtom('CN' , 'CT2',  -0.18 ,  10,  -1),
   31: CharmmAtom('HN1', 'HA' ,   0.09 ,  10,  -1),
   32: CharmmAtom('HN2', 'HA' ,   0.09 ,  10,  -1),   # end group
    4: CharmmAtom('NG' , 'NY' ,  -0.03 ,  11,  -1),   # C->N
    3: CharmmAtom('CD1', 'CA' ,   0.035,  11,  -1),
   30: CharmmAtom('HD1', 'HP' ,   0.115,  11,  -1),
    2: CharmmAtom('NE1', 'NY' ,  -0.61 ,  11,  -1),
    1: CharmmAtom('CX1', 'CT3',   0.11 ,  11,  -1),   # added methyl
   27: CharmmAtom('HX1', 'HA',    0.09 ,  11,  -1),   # added methyl
   28: CharmmAtom('HX2', 'HA',    0.09 ,  11,  -1),   # added methyl
   29: CharmmAtom('HX3', 'HA',    0.09 ,  11,  -1),   # added methyl
   21: CharmmAtom('CE2', 'CPT',   0.13 ,  11,  -1),
   22: CharmmAtom('CD2', 'CPT',  -0.02 ,  11,  -1),   # end group
   23: CharmmAtom('CE3', 'CA' ,  -0.115,  12,  -1),
   51: CharmmAtom('HE3', 'HP' ,   0.115,  12,  -1),   # end group
   24: CharmmAtom('CZ3', 'CA' ,  -0.115,  13,  -1),
   52: CharmmAtom('HZ3', 'HP' ,   0.115,  13,  -1),   # end group
   25: CharmmAtom('CZ2', 'CA' ,  -0.115,  14,  -1),
   53: CharmmAtom('HZ2', 'HP' ,   0.115,  14,  -1),   # end group
   26: CharmmAtom('CH2', 'CA' ,  -0.115,  15,  -1),
   54: CharmmAtom('HH2', 'HP' ,   0.115,  15,  -1),   # end group
}

mol = Chem.MolFromSmiles('CN1C=[N+](CCCCCCOC(=O)CC[C@H](N)C(O)=O)C2=C1C=CC=C2')
mol_h = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol_h)
AllChem.MMFFOptimizeMolecule(mol_h)
Chem.MolToPDBFile(mol_h, 'phlg_monomer.pdb')

group_map = defaultdict(list)
for atom_serial, atom_data in atom_map.items():
    new_atom_data = CharmmAtom(
        atom_data.label,
        atom_data.ff_type,
        atom_data.charge,
        atom_data.charge_group_id,
        atom_serial
    )
    group_map[atom_data.charge_group_id].append(new_atom_data)

with open('top_ligand.inp', 'w') as f:
    for group_id in sorted(list(group_map.keys())):
        group_charge = sum(atom.charge for atom in group_map[group_id])
        f.write('GROUP  ! %.2f\n' % group_charge)
        atom_list = group_map[group_id]
        sorted_atom_list = sorted(atom_list, key=lambda a: a.serial)
        for atom in sorted_atom_list:
            f.write('ATOM  %3s   %3s   %5.2f\n' % (
                atom.label,
                atom.ff_type,
                atom.charge
            ))

    row = 'BOND'
    for bond in mol_h.GetBonds():
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        begin_atom_serial = begin_atom.GetIdx() + 1
        end_atom_serial = end_atom.GetIdx() + 1
        row += '    %-3s  %-3s' % (
            atom_map[begin_atom_serial].label,
            atom_map[end_atom_serial].label
        )
        if (bond.GetIdx() + 1) % 4 == 0:
            f.write(row + '\n')
            row = 'BOND'
    if row != 'BOND':
        f.write(row + '\n')

# Use this to replace atom labels in PDB file
with open('atom_list.txt', 'w') as f:
    for atom_serial in sorted(list(atom_map.keys())):
        f.write('%-3d  %-4s\n' % (
            atom_serial,
            atom_map[atom_serial].label
        ))

edge_list = []
for bond in mol_h.GetBonds():
    begin_atom = bond.GetBeginAtom()
    end_atom = bond.GetEndAtom()
    edge_list.append((
        begin_atom.GetIdx(),
        end_atom.GetIdx()
    ))
graph = nx.Graph(edge_list)
for atom in mol_h.GetAtoms():
    node_key = atom.GetIdx()
    symbol = atom.GetSymbol()
    graph.nodes[node_key]['symbol'] = symbol
with open('molecule.json', 'w') as f:
    json.dump(nx.node_link_data(graph), f, indent=2)
