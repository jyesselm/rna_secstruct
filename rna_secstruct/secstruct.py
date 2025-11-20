"""
representation of secondary structure with motif
"""

from typing import List, Optional, Union
from dataclasses import dataclass
from rna_secstruct.parser import Parser, is_valid_dot_bracket_str
from rna_secstruct.motif import Motif


@dataclass(order=True)
class MotifSearchParams:
    """
    params for rna design
    """

    sequence: Optional[str] = None
    structure: Optional[str] = None
    m_type: Optional[str] = None
    min_pos: int = 0
    max_pos: int = 999
    min_id: int = 0
    max_id: int = 999


class SecStruct:
    """
    A class to represent a secondary structure.
    """

    def __init__(self, sequence: str, structure: str):
        """Initialize a SecStruct with sequence and structure.

        Parsing is deferred until motifs are actually accessed (lazy loading).
        This improves performance for large structures when only basic operations
        are needed.

        Args:
            sequence: A sequence of nucleotides.
            structure: A dot bracket structure.

        Raises:
            ValueError: If sequence and structure have different lengths or different number of strands.
        """
        # Validate inputs but don't parse yet
        if len(sequence) != len(structure):
            raise ValueError("sequence and structure must be the same length")
        if sequence.count("&") != structure.count("&"):
            raise ValueError("sequence and structure must have the same number of strands")

        self.__sequence = sequence
        self.__structure = structure
        # Lazy loading: parse only when motifs are accessed
        self.__root = None
        self.__motifs = None
        self.__parsed = False

    def __len__(self) -> int:
        """Return the length of the sequence.

        Returns:
            int: Length of sequence/structure.
        """
        return len(self.__sequence)

    def __add__(self, other):
        """Add two secondary structures together.

        Args:
            other: Another SecStruct to add.

        Returns:
            SecStruct: A new SecStruct instance with concatenated sequence and structure.
        """
        return SecStruct(self.__sequence + other.__sequence, self.__structure + other.__structure)

    def itermotifs(self):
        """Iterate over motifs.

        Returns:
            Iterator: Iterator over (motif_id, motif) tuples.
        """
        return iter(self.motifs.items())

    def __iter__(self):
        """Iterate over the motifs.

        Returns:
            Iterator: Iterator over motif values.
        """
        return iter(self.motifs.values())

    def __getitem__(self, item) -> Union['SecStruct', Motif]:
        """Support both slicing (returns new SecStruct) and motif access (returns Motif).

        Slicing follows immutable container pattern - returns new instance.

        Args:
            item: Slice object or motif ID.

        Returns:
            Union[SecStruct, Motif]: New SecStruct for slicing, Motif for ID access.

        Raises:
            ValueError: If motif ID is not found.
        """
        if isinstance(item, slice):
            # Slicing returns new SecStruct (immutable)
            return SecStruct(
                self.__sequence[item],
                self.__structure[item]
            )
        # Existing motif access logic (returns Motif, not SecStruct)
        if item not in self.motifs:
            raise ValueError(f"no motif with id {item}")
        return self.motifs[item]

    def __repr__(self):
        """String representation - does not trigger parsing"""
        return f"{self.__sequence}, {self.__structure}"

    def __parse(self):
        """Internal: Parse structure when needed (lazy loading).

        This is called automatically when motifs are first accessed.
        """
        if not self.__parsed:
            self.__root = Parser().parse(self.__sequence, self.__structure)
            self.__motifs = self.__get_motifs(self.__root)
            self.__parsed = True

    def __get_motifs(self, root):
        """Internal: Get a dictionary of motifs from root.

        Args:
            root: Root motif to collect from.

        Returns:
            Dict[int, Motif]: Dictionary mapping motif IDs to Motif objects.
        """
        motifs = {}
        self.__collect_motifs(root, motifs)
        return motifs

    def __collect_motifs(self, cur_motif, motifs):
        """Internal: Recursively collect motifs from motif tree.

        Args:
            cur_motif: Current motif to process.
            motifs: Dictionary to populate with collected motifs.
        """
        motifs[cur_motif.m_id] = cur_motif
        for c in cur_motif.children:
            self.__collect_motifs(c, motifs)

    # design ###################################################################

    def change_motif_sequence(self, m_id, seqeuence):
        """Change the sequence of a motif.

        Args:
            m_id: The motif ID.
            seqeuence: The new sequence.

        Note:
            This method is not yet implemented.
        """
        pass

    def change_motif(self, m_id, sequence, structure) -> None:
        """Change a motif sequence and secondary structure triggering a reparse.

        NOTE: This method modifies the SecStruct in place. For immutable
        operations, consider using replace_motif() which returns a new instance.

        Args:
            m_id: The ID of the motif to change.
            sequence: The new sequence.
            structure: The new structure.

        Raises:
            ValueError: If sequence and structure have different lengths,
                different number of strands, or if strands cannot be added to motif.
        """
        # Ensure we're parsed before modifying
        if not self.__parsed:
            self.__parse()

        is_valid_dot_bracket_str(structure)
        if len(sequence) != len(structure):
            raise ValueError("sequence and structure must be the same length")
        if sequence.count("&") != structure.count("&"):
            raise ValueError("sequence and structure must have the same number of strands")
        m = self[m_id]
        # logic to change motif type only can in number of strands
        if m.num_strands() < sequence.count("&") + 1:
            raise ValueError("cannot add strands to a motif")
        m_type = m.m_type
        if m.num_strands() > sequence.count("&") + 1 and sequence.count("&") == 0:
            if structure.count("(") == 0:
                m_type = "SINGLESTRAND"
            else:
                m_type = "HAIRPIN"
        if m.has_parent() and m.parent.is_helix():
            self.__change_inner_flanking(m.parent, sequence[0] + sequence[-1])
        strands = sequence.split("&")
        for s1, s2, c in zip(strands[:-1], strands[1:], m.children):
            self.__change_outer_flanking(c, s1[-1] + s2[0])
        m.sequence = sequence
        m.structure = structure
        m.m_type = m_type
        full_seq = self.__root.recursive_sequence()
        full_ss = self.__root.recursive_structure()
        self.__sequence = full_seq
        self.__structure = full_ss
        self.__root = Parser().parse(full_seq, full_ss)
        self.__motifs = self.__get_motifs(self.__root)

    def __change_inner_flanking(self, m, new_cp: str):
        """Internal: Change the inner flanking sequence of a motif.

        Args:
            m: The motif to modify.
            new_cp: New complementary pair string (2 characters).
        """
        it = m.sequence.find("&")
        tks = list(m.sequence)
        tks[it - 1], tks[it + 1] = new_cp[0], new_cp[1]
        m.sequence = "".join(tks)

    def __change_outer_flanking(self, m, new_cp: str):
        """Internal: Change the outer flanking sequence of a motif.

        Args:
            m: The motif to modify.
            new_cp: New complementary pair string (2 characters).
        """
        tks = list(m.sequence)
        tks[0], tks[-1] = new_cp[0], new_cp[1]
        m.sequence = "".join(tks)

    # search ###################################################################

    def __get_motifs_by_params(self, msp: MotifSearchParams) -> List[Motif]:
        """Internal: Get a list of motifs matching search parameters.

        Args:
            msp: Motif search parameters.

        Returns:
            List[Motif]: List of motifs matching the search criteria.
        """
        motifs = []
        for m in self:
            if msp.m_type is not None and m.m_type != msp.m_type:
                continue
            if msp.sequence is not None and msp.sequence != m.sequence:
                continue
            if msp.structure is not None and msp.structure != m.structure:
                continue
            if m.m_id < msp.min_id or m.m_id > msp.max_id:
                continue
            if m.start_pos < msp.min_pos or m.end_pos > msp.max_pos:
                continue
            motifs.append(m)
        return motifs

    def get_motifs_by_strand_lengths(
        self, strand_lengths, msp: Optional[MotifSearchParams] = None
    ) -> List[Motif]:
        """Get motifs by the length of each strand.

        Strand lengths is a list of n size, where n is the size of each strand.
        The list is in order of the strand size. For example, if the motif is a
        hairpin it can only have 1 strand, so the list will be [n]. If the motif
        is a junction it can have 2 or more strands, so the list will be [n, m, ...].

        Args:
            strand_lengths: The list of strand lengths, e.g. [3, 4, 5].
            msp: Optional motif search params.

        Returns:
            List[Motif]: List of motifs matching the strand lengths.
        """
        if msp is None:
            msp = MotifSearchParams()
        selected_motifs = self.__get_motifs_by_params(msp)
        motifs = []
        for m in selected_motifs:
            # get the length of each strand
            lengths = [len(s) for s in m.strands]
            # compare the lengths
            if lengths == strand_lengths:
                motifs.append(m)
        return motifs

    def get_twoway_junctions_by_topology(
        self, x_pos, y_pos, msp: Optional[MotifSearchParams] = None
    ) -> List[Motif]:
        """Get two-way junctions by topology.

        Args:
            x_pos: X position in topology.
            y_pos: Y position in topology.
            msp: Optional motif search params. If provided, m_type must be JUNCTION or None.

        Returns:
            List[Motif]: List of two-way junction motifs matching the topology.

        Raises:
            ValueError: If msp.m_type is provided and is not JUNCTION.
        """
        # add 2 two each number in topology
        # to account for the 2 flanking base pairs
        if msp is None:
            msp = MotifSearchParams(m_type="JUNCTION")
        else:
            if msp.m_type is not None and msp.m_type != "JUNCTION":
                raise ValueError("m_type must be JUNCTION")
            msp.m_type = "JUNCTION"

        topology = [t + 2 for t in [x_pos, y_pos]]
        return self.get_motifs_by_strand_lengths(topology, msp)

    def get_motifs(self, msp) -> List[Motif]:
        """Get motifs by sequence and structure.

        Args:
            msp: Motif search parameters.

        Returns:
            List[Motif]: List of motifs matching the search criteria.
        """
        return self.__get_motifs_by_params(msp)

    def get_motifs_by_token(self, token, msp: Optional[MotifSearchParams] = None) -> List[Motif]:
        """Get a list of motifs by a token.

        Args:
            token: The token to search for.
            msp: Optional motif search params.

        Returns:
            List[Motif]: List of motifs with matching token.
        """
        if msp is None:
            msp = MotifSearchParams()
        selected_motifs = self.__get_motifs_by_params(msp)
        motifs = []
        for m in selected_motifs:
            if token == m.token:
                motifs.append(m)
        return motifs

    # properites ###############################################################
    @property
    def _root(self):
        """Get the root motif (lazy-loaded).

        This property triggers parsing if not already done.

        Returns:
            Motif: The root motif of the structure.
        """
        if not self.__parsed:
            self.__parse()
        return self.__root

    @property
    def motifs(self):
        """Get motifs dictionary (lazy-loaded).

        This property triggers parsing if not already done.

        Returns:
            Dict[int, Motif]: Dictionary mapping motif IDs to Motif objects.
        """
        if not self.__parsed:
            self.__parse()
        return self.__motifs

    @property
    def sequence(self):
        """Get the sequence.

        Returns:
            str: The RNA sequence.
        """
        return self.__sequence

    @property
    def structure(self):
        """Get the structure.

        Returns:
            str: The RNA secondary structure.
        """
        return self.__structure

    # getters and setters ######################################################
    def get_copy(self):
        """Get a copy of this secondary structure.

        Returns:
            SecStruct: A new SecStruct instance with the same sequence and structure.
        """
        return SecStruct(self.__sequence, self.__structure)

    def get_hairpins(self):
        """Get hairpin motifs.

        Returns:
            List[Motif]: List of hairpin motifs.
        """
        return [m for m in self if m.is_hairpin()]

    def get_helices(self):
        """Get helix motifs.

        Returns:
            List[Motif]: List of helix motifs.
        """
        return [m for m in self if m.is_helix()]

    def get_junctions(self):
        """Get junction motifs.

        Returns:
            List[Motif]: List of junction motifs.
        """
        return [m for m in self if m.is_junction()]

    def get_single_strands(self):
        """Get single strand motifs.

        Returns:
            List[Motif]: List of single strand motifs.
        """
        return [m for m in self if m.is_single_strand()]

    def get_num_motifs(self):
        """Get the number of motifs.

        Returns:
            int: The number of motifs in the structure.
        """
        # Access motifs property to trigger lazy parsing if needed
        return len(self.motifs)

    def get_sub_structure(self, root_id):
        """Get a substructure starting from a given motif.

        Args:
            root_id: The motif ID to use as the root of the substructure.

        Returns:
            SecStruct: A new SecStruct instance representing the substructure.
        """
        m = self[root_id]
        seq = m.recursive_sequence()
        struct = m.recursive_structure()
        # fix dropping the first basepair in the structure
        if m.has_parent() and m.parent.is_helix():
            seq_spl = m.parent.sequence.split("&")
            seq = seq_spl[0][-1] + seq + seq_spl[1][0]
            struct = "(" + struct + ")"
        return SecStruct(seq, struct)

    def to_str(self):
        """Get a string representation of this secondary structure.

        Returns:
            str: String representation of the structure.
        """
        # Access _root property to trigger lazy parsing if needed
        return self._root.to_str()

    # SequenceStructure-like methods (immutable container pattern) ############

    def split_strands(self) -> List['SecStruct']:
        """Split both sequence and structure over '&' and return list of SecStruct objects.

        Returns new SecStruct instances (immutable container pattern).

        Returns:
            List[SecStruct]: List of new SecStruct objects, one per strand.
        """
        seqs = self.__sequence.split("&")
        structs = self.__structure.split("&")
        return [SecStruct(s, st) for s, st in zip(seqs, structs)]

    def insert(self, pos: int, other: 'SecStruct') -> 'SecStruct':
        """Insert a SecStruct object at a given position.

        Returns a NEW SecStruct instance. Original is unchanged (immutable container).

        Args:
            pos: Position to insert at.
            other: SecStruct to insert.

        Returns:
            SecStruct: New SecStruct instance with insertion.

        Raises:
            ValueError: If position is invalid.
        """
        if pos < 0 or pos > len(self.__sequence):
            raise ValueError(f"Invalid position: {pos}. Must be between 0 and {len(self.__sequence)}")
        seq = self.__sequence[:pos] + other.__sequence + self.__sequence[pos:]
        struct = self.__structure[:pos] + other.__structure + self.__structure[pos:]
        return SecStruct(seq, struct)

    def join(self, other: 'SecStruct') -> 'SecStruct':
        """Join two SecStruct objects with '&' separating each strand.

        Returns a NEW SecStruct instance. Original is unchanged (immutable container).

        Args:
            other: SecStruct to join.

        Returns:
            SecStruct: New SecStruct instance.
        """
        return SecStruct(
            self.__sequence + "&" + other.__sequence,
            self.__structure + "&" + other.__structure,
        )

    def replace(self, other: 'SecStruct', pos: int) -> 'SecStruct':
        """Replace sequence and structure at specified position.

        Returns a NEW SecStruct instance. Original is unchanged (immutable container).

        Args:
            other: SecStruct to replace with.
            pos: Position to replace at.

        Returns:
            SecStruct: New SecStruct instance with replacement.

        Raises:
            ValueError: If position is invalid or replacement extends beyond structure length.
        """
        if pos < 0 or pos > len(self.__sequence):
            raise ValueError(f"Invalid position: {pos}. Must be between 0 and {len(self.__sequence)}")
        if pos + len(other.__sequence) > len(self.__sequence):
            raise ValueError(
                f"Replacement extends beyond structure length. "
                f"Position {pos} + length {len(other.__sequence)} > structure length {len(self.__sequence)}"
            )
        sequence = (
            self.__sequence[:pos]
            + other.__sequence
            + self.__sequence[pos + len(other.__sequence):]
        )
        structure = (
            self.__structure[:pos]
            + other.__structure
            + self.__structure[pos + len(other.__structure):]
        )
        return SecStruct(sequence, structure)

    def remove(self, start: int, end: int) -> 'SecStruct':
        """Remove a region from the structure.

        Returns a NEW SecStruct instance. Original is unchanged (immutable container).

        Args:
            start: Start position (inclusive).
            end: End position (exclusive).

        Returns:
            SecStruct: New SecStruct instance with region removed.

        Raises:
            ValueError: If range is invalid.
        """
        if start < 0 or end > len(self.__sequence) or start >= end:
            raise ValueError(
                f"Invalid range: {start} to {end}. "
                f"Start must be >= 0, end must be <= {len(self.__sequence)}, and start < end"
            )
        seq = self.__sequence[:start] + self.__sequence[end:]
        struct = self.__structure[:start] + self.__structure[end:]
        return SecStruct(seq, struct)

    def subtract(self, other: 'SecStruct') -> 'SecStruct':
        """Remove a substructure from this structure (if found).

        Returns a NEW SecStruct instance. Original is unchanged (immutable container).

        Args:
            other: SecStruct to remove.

        Returns:
            SecStruct: New SecStruct instance with substructure removed.

        Raises:
            ValueError: If substructure is not found.
        """
        # Find position of other in self
        pos = self.__sequence.find(other.__sequence)
        if pos == -1:
            raise ValueError(
                f"Substructure not found. "
                f"Sequence '{other.__sequence}' not found in '{self.__sequence}'"
            )
        return self.remove(pos, pos + len(other.__sequence))

    def to_dict(self) -> dict:
        """Return a dictionary representation of the SecStruct.

        Returns:
            dict: Dictionary with 'sequence' and 'structure' keys.
        """
        return {"sequence": self.__sequence, "structure": self.__structure}

    def to_comma_delimited(self) -> str:
        """Return a CSV representation of the SecStruct.

        Returns:
            str: Comma-delimited string: "sequence,structure".
        """
        return f"{self.__sequence},{self.__structure}"
