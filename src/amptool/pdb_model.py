from collections import namedtuple
from collections.abc import Callable
from pathlib import Path

from amptool.models import Protein, ProteinResidue, ProteinChain, ProteinAtom, IntCoor

AtomRow = namedtuple(
    "AtomRow",
    [
        "record_label",
        "atom_serial",
        "atom_label",
        "alternate_location",
        "residue_label",
        "chain_id",
        "residue_serial",
        "insertion_code",
        "coordinate_x",
        "coordinate_y",
        "coordinate_z",
        "occupancy",
        "temperature_factor",
        "segment_id",
        "element_symbol",
        "atom_charge",
    ],
)


def make_row(row: str) -> AtomRow | None:
    """Parse an ATOM or HETATM PDB record into an AtomRow; return None for other records."""
    if row.startswith("ATOM") or row.startswith("HETATM"):
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
            make_field(row, 78, 80, float),
        )


def make_field(
    row: str, begin: int, end: int, converter: Callable[[str], str | int | float]
) -> str | int | float | None:
    """Extract a fixed-width field from a PDB record and convert it; return None on failure."""
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


class ProteinPDBBuilder(object):
    """Parses a PDB file into a Protein object."""

    def __init__(self, pdb_file: str) -> None:
        self.__pdb_file = pdb_file

    def build(self) -> Protein:
        """Read the PDB file and return a fully assembled Protein."""
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
                        chain_residues.append(
                            ProteinResidue(
                                previous_residue_serial,
                                previous_residue_label,
                                residue_atoms,
                            )
                        )
                        residue_atoms = []

                if atom_row.chain_id != previous_chain_id:
                    if chain_residues:
                        protein_chains.append(
                            ProteinChain(previous_chain_id, chain_residues)
                        )
                        chain_residues = []
                residue_atoms.append(
                    ProteinAtom(
                        atom_row.atom_serial,
                        atom_row.atom_label,
                        atom_row.coordinate_x,
                        atom_row.coordinate_y,
                        atom_row.coordinate_z,
                    )
                )
                previous_chain_id = atom_row.chain_id
                previous_residue_label = atom_row.residue_label
                previous_residue_serial = atom_row.residue_serial

        if residue_atoms:
            chain_residues.append(
                ProteinResidue(
                    previous_residue_serial, previous_residue_label, residue_atoms
                )
            )

        if chain_residues:
            protein_chains.append(ProteinChain(previous_chain_id, chain_residues))

        return Protein(protein_chains)


def make_protein_model(monomer_ic_path: Path, protein_path: Path) -> None:
    """

    :param monomer_ic_path:
    :param protein_path:
    :return:
    """
    with monomer_ic_path.open("r") as f:
        ic_table = []
        for row in f.readlines():
            ic_table.append(IntCoor(*row.split()))

    chain_path = Path(__file__).parent / "assets/all_ala_helix_short.pdb"
    chain_str = chain_path.resolve().as_posix()

    protein = ProteinPDBBuilder(chain_str).build()
    protein.remove_side_chains()
    protein.fill_side_chains(ic_table)
    protein.renumerate_atoms()

    with protein_path.open("w") as f:
        protein.to_pdb(f)
