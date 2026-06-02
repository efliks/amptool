import math
from collections import OrderedDict
from typing import IO

import numpy


class ProteinAtom(object):
    """A single atom with a serial number, label, and Cartesian coordinates."""

    def __init__(
            self,
            atom_serial: int,
            atom_label: str,
            location_x: float,
            location_y: float,
            location_z: float,
    ) -> None:
        self.serial = atom_serial
        self.label = atom_label
        self.x = location_x
        self.y = location_y
        self.z = location_z


class IntCoor(object):
    """Internal coordinate entry: atom labels and the bond length, angle, and dihedral
    used to reconstruct the position of atom i from atoms j, k, and l."""

    def __init__(
            self,
            i: str,
            j: str,
            k: str,
            l: str,
            ij: str | float,
            ijk: str | float,
            ijkl: str | float,
    ) -> None:
        self.i = i
        self.j = j
        self.k = k
        self.l = l
        self.ij = float(ij)
        self.ijk = float(ijk)
        self.ijkl = float(ijkl)

    @staticmethod
    def __normalize_vector(vector: numpy.ndarray) -> numpy.ndarray:
        """Return a unit vector in the direction of vector."""
        return (1.0 / numpy.linalg.norm(vector)) * vector

    @staticmethod
    def __make_vector(atom: ProteinAtom) -> numpy.ndarray:
        """Return the Cartesian position of atom as a numpy array."""
        return numpy.array((atom.x, atom.y, atom.z))

    def make_atom(
            self, atom_map: dict[str, ProteinAtom], atom_serial: int = 0
    ) -> ProteinAtom:
        """Reconstruct atom i's 3D position from internal coordinates and atoms j, k, l."""
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
        return ProteinAtom(atom_serial, self.i, location[0], location[1], location[2])


class ProteinResidue(object):
    """An amino acid residue containing an ordered map of its atoms."""

    BACKBONE_LABELS = ("N", "HN", "CA", "HA", "C", "O", "CB")

    def __init__(
            self,
            residue_serial: int,
            residue_label: str,
            atom_list: list[ProteinAtom],
    ) -> None:
        self.serial = residue_serial
        self.label = residue_label
        self.atom_map = OrderedDict(
            [(protein_atom.label, protein_atom) for protein_atom in atom_list]
        )

    def __iter__(self):
        for atom in self.atom_map.values():
            yield atom

    def remove_side_chain(self) -> None:
        """Discard all atoms except the backbone labels."""
        self.atom_map = OrderedDict(
            [
                (protein_atom.label, protein_atom)
                for protein_atom in self.atom_map.values()
                if protein_atom.label in ProteinResidue.BACKBONE_LABELS
            ]
        )

    def fill_side_chain(self, ic_table: list[IntCoor]) -> None:
        """Append side-chain atoms computed from internal coordinates."""
        for internal_coordinate in ic_table:
            atom = internal_coordinate.make_atom(self.atom_map)
            self.atom_map[atom.label] = atom

    def renumerate_atoms(self, atom_serial: int) -> int:
        """Reassign serial numbers sequentially from atom_serial; return the next serial."""
        reordered_map = []
        for atom in self:
            reordered_map.append(
                (
                    atom.label,
                    ProteinAtom(atom_serial, atom.label, atom.x, atom.y, atom.z),
                )
            )
            atom_serial += 1
        self.atom_map = OrderedDict(reordered_map)
        return atom_serial


class ProteinChain(object):
    """A polypeptide chain containing an ordered map of its residues."""

    def __init__(self, chain_id: str, residue_list: list[ProteinResidue]) -> None:
        self.chain_id = chain_id
        self.residue_map = OrderedDict(
            [
                (protein_residue.serial, protein_residue)
                for protein_residue in residue_list
            ]
        )

    def __iter__(self):
        for residue in self.residue_map.values():
            yield residue


class Protein(object):
    """A complete protein containing an ordered map of its chains."""

    def __init__(self, chain_list: list[ProteinChain]) -> None:
        self.chain_map = OrderedDict(
            [(protein_chain.chain_id, protein_chain) for protein_chain in chain_list]
        )

    def __iter__(self):
        for chain in self.chain_map.values():
            yield chain

    def remove_side_chains(self) -> None:
        """Remove side-chain atoms from every residue in the protein."""
        for chain in self:
            for residue in chain:
                residue.remove_side_chain()

    def fill_side_chains(self, ic_table: list[IntCoor]) -> None:
        """Fill side-chain atoms for every residue using internal coordinates."""
        for chain in self:
            for residue in chain:
                residue.fill_side_chain(ic_table)

    def renumerate_atoms(self) -> None:
        """Reassign all atom serial numbers sequentially across every chain and residue."""
        atom_serial = 1
        for chain in self:
            for residue in chain:
                atom_serial = residue.renumerate_atoms(atom_serial)

    def to_pdb(self, fp: IO[str]) -> None:
        """Write ATOM records for all atoms followed by an END record to fp."""
        for chain in self:
            for residue in chain:
                for atom in residue:
                    fp.write(
                        f"ATOM  {atom.serial:>5}  {atom.label:<4}"
                        f"{residue.label}{chain.chain_id}   {residue.serial:>3}    "
                        f"{atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}\n"
                    )
        fp.write("END\n")
