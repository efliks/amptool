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
            return converter(column)
        except ValueError:
            pass


class ProteinAtom(object):
    def __init__(self, atom_serial, atom_label,
                 location_x, location_y, location_z):
        self.atom_serial = atom_serial
        self.atom_name = atom_label
        self.location_x = location_x
        self.location_y = location_y
        self.location_z = location_z


class ProteinResidue(object):
    def __init__(self, residue_serial, residue_label, atom_list):
        self.residue_serial = residue_serial
        self.residue_name = residue_label
        self.atom_map = OrderedDict([
            (protein_atom.atom_name, protein_atom)
            for protein_atom in atom_list
        ])


class ProteinChain(object):
    def __init__(self, chain_id, residue_list):
        self.chain_id = chain_id
        self.residue_map = OrderedDict([
            (protein_residue.residue_serial, protein_residue)
            for protein_residue in residue_list
        ])


class Protein(object):
    def __init__(self, chain_list):
        self.chain_map = OrderedDict([
            (protein_chain.chain_id, protein_chain)
            for protein_chain in chain_list
        ])


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


protein = ProteinPDBBuilder('all_ala_helix.pdb').build()
print('!')
