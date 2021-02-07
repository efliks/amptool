import math
import numpy

from collections import namedtuple, OrderedDict


AtomRow = namedtuple('AtomRow', [
    'record_label',
    'atom_serial',
    'atom_label',
    'alternate_location',
    'residue_label',
    'chain_id',
    'residue_serial',
    'insertion_code',
    'coordinate_x',
    'coordinate_y',
    'coordinate_z',
    'occupancy',
    'temperature_factor',
    'segment_id',
    'element_symbol',
    'atom_charge'
])


def make_row(row):
    if row.startswith('ATOM') or row.startswith('HETATM'):
        return AtomRow(
            make_field(row, 0, 6, str),
            make_field(row, 6, 11, int),
            make_field(row, 12, 16, str),
            make_field(row, 16, 17, str),
            make_field(row, 17, 20, str),
            make_field(row, 21, -1, str),
            make_field(row, 22, 26, int),
            make_field(row, 26, -1, str),
            make_field(row, 30, 38, float),
            make_field(row, 38, 46, float),
            make_field(row, 46, 54, float),
            make_field(row, 54, 60, float),
            make_field(row, 60, 66, float),
            make_field(row, 72, 75, str),
            make_field(row, 76, 78, str),
            make_field(row, 78, 80, float)
        )


def make_field(row, begin, end, converter):
    if begin < end:
        column = row[begin:end]
    else:
        column = row[begin]
    if column:
        try:
            converted = converter(column)
            if isinstance(converted, str):
                return converted.strip()
            return converted
        except ValueError:
            pass


class ProteinAtom(object):
    def __init__(self, atom_serial, atom_label,
                 location_x, location_y, location_z):
        self.serial = atom_serial
        self.label = atom_label
        self.x = location_x
        self.y = location_y
        self.z = location_z


class IntCoor(object):
    def __init__(self, i, j, k, l, ij, ijk, ijkl):
        self.i = i
        self.j = j
        self.k = k
        self.l = l
        self.ij = float(ij)
        self.ijk = float(ijk)
        self.ijkl = float(ijkl)

    @staticmethod
    def __normalize_vector(vector):
        return (1.0 / numpy.linalg.norm(vector)) * vector

    @staticmethod
    def __make_vector(atom):
        return numpy.array((
            atom.x,
            atom.y,
            atom.z
        ))

    def make_atom(self, atom_map, atom_serial=0):
        #         (a) I
        #              \
        #               \
        #            (b) J----K (c)
        #                      \
        #                       \
        #                        L (d)
        #     values (Rij),(Tijk),(Pijkl),(Tjkl),(Rkl)
        a = self.__make_vector(atom_map[self.l])
        b = self.__make_vector(atom_map[self.k])
        c = self.__make_vector(atom_map[self.j])
        ba = a - b
        j = self.__normalize_vector(c - b)
        k = self.__normalize_vector(numpy.cross(ba, j))
        i = self.__normalize_vector(numpy.cross(j, k))
        psi = self.ijkl * math.pi / 180.0
        t = math.cos(psi) * i + math.sin(psi) * k
        chi = self.ijk * math.pi / 180.0
        q = -math.cos(chi) * j + math.sin(chi) * t
        location = c + self.ij * q
        return ProteinAtom(
            atom_serial,
            self.i,
            location[0],
            location[1],
            location[2]
        )


class ProteinResidue(object):
    BACKBONE_LABELS = ('N', 'HN', 'CA', 'HA', 'C', 'O', 'CB')

    def __init__(self, residue_serial, residue_label, atom_list):
        self.serial = residue_serial
        self.label = residue_label
        self.atom_map = OrderedDict([
            (protein_atom.label, protein_atom)
            for protein_atom in atom_list
        ])

    def __iter__(self):
        for atom in self.atom_map.values():
            yield atom

    def remove_side_chain(self):
        self.atom_map = OrderedDict([
            (protein_atom.label, protein_atom)
            for protein_atom in self.atom_map.values()
            if protein_atom.label in ProteinResidue.BACKBONE_LABELS
        ])

    def fill_side_chain(self, ic_table):
        for internal_coordinate in ic_table:
            atom = internal_coordinate.make_atom(self.atom_map)
            self.atom_map[atom.label] = atom

    def renumerate_atoms(self, atom_serial):
        reordered_map = []
        for atom in self:
            reordered_map.append((atom.label, ProteinAtom(
                atom_serial,
                atom.label,
                atom.x,
                atom.y,
                atom.z
            )))
            atom_serial += 1
        self.atom_map = OrderedDict(reordered_map)
        return atom_serial


class ProteinChain(object):
    def __init__(self, chain_id, residue_list):
        self.chain_id = chain_id
        self.residue_map = OrderedDict([
            (protein_residue.serial, protein_residue)
            for protein_residue in residue_list
        ])

    def __iter__(self):
        for residue in self.residue_map.values():
            yield residue


class Protein(object):
    def __init__(self, chain_list):
        self.chain_map = OrderedDict([
            (protein_chain.chain_id, protein_chain)
            for protein_chain in chain_list
        ])

    def __iter__(self):
        for chain in self.chain_map.values():
            yield chain

    def remove_side_chains(self):
        for chain in self:
            for residue in chain:
                residue.remove_side_chain()

    def fill_side_chains(self, ic_table):
        for chain in self:
            for residue in chain:
                residue.fill_side_chain(ic_table)

    def renumerate_atoms(self):
        atom_serial = 1
        for chain in self:
            for residue in chain:
                atom_serial = residue.renumerate_atoms(
                    atom_serial)

    def to_pdb(self, fp):
        for chain in self:
            for residue in chain:
                for atom in residue:
                    fp.write(
                        f'ATOM  {atom.serial:>5}  {atom.label:<4}'
                        f'{residue.label}{chain.chain_id}   {residue.serial:>3}    '
                        f'{atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}\n'
                    )
        fp.write('END\n')


class ProteinPDBBuilder(object):
    def __init__(self, pdb_file):
        self.__pdb_file = pdb_file

    def build(self):
        with open(self.__pdb_file) as f:
            previous_chain_id = None
            previous_residue_label = None
            previous_residue_serial = None
            protein_chains = []
            chain_residues = []
            residue_atoms = []
            for row in f.readlines():
                atom_row = make_row(row)
                if not atom_row:
                    continue
                if atom_row.residue_serial != previous_residue_serial:
                    if residue_atoms:
                        chain_residues.append(ProteinResidue(
                            previous_residue_serial,
                            previous_residue_label,
                            residue_atoms
                        ))
                        residue_atoms = []
                if atom_row.chain_id != previous_chain_id:
                    if chain_residues:
                        protein_chains.append(ProteinChain(
                            previous_chain_id,
                            chain_residues
                        ))
                        chain_residues = []
                residue_atoms.append(ProteinAtom(
                    atom_row.atom_serial,
                    atom_row.atom_label,
                    atom_row.coordinate_x,
                    atom_row.coordinate_y,
                    atom_row.coordinate_z
                ))
                previous_chain_id = atom_row.chain_id
                previous_residue_label = atom_row.residue_label
                previous_residue_serial = atom_row.residue_serial
        if residue_atoms:
            chain_residues.append(ProteinResidue(
                previous_residue_serial,
                previous_residue_label,
                residue_atoms
            ))
        if chain_residues:
            protein_chains.append(ProteinChain(
                previous_chain_id,
                chain_residues
            ))
        return Protein(protein_chains)


with open('monomer_reorder.ic') as f:
    ic_table = []
    for row in f.readlines():
        ic_table.append(IntCoor(*row.split()))

protein = ProteinPDBBuilder('all_ala_helix_short.pdb').build()
protein.remove_side_chains()
protein.fill_side_chains(ic_table)
protein.renumerate_atoms()

with open('my.pdb', 'w') as f:
    protein.to_pdb(f)
