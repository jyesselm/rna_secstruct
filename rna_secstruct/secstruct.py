"""
representation of secondary structure with motif
"""

from typing import List, Optional, Union, Tuple
from dataclasses import dataclass
from rna_secstruct.parser import Parser, is_valid_dot_bracket_str
from rna_secstruct.motif import Motif
from rna_secstruct.connectivity import get_connectivity_list


@dataclass(order=True)
class MotifSearchParams:
    """Parameters for RNA motif search.

    Attributes:
        sequence: Exact sequence to match.
        structure: Exact structure to match.
        m_type: Motif type to match (e.g., 'HELIX', 'HAIRPIN', 'JUNCTION').
        min_pos: Minimum start position.
        max_pos: Maximum end position.
        min_id: Minimum motif ID.
        max_id: Maximum motif ID.
        token: Token to match.
        min_length: Minimum motif length.
        max_length: Maximum motif length.
        strand_lengths: List of strand lengths to match.
        has_children: Whether motif must have children (True) or not (False), or None for any.
    """

    sequence: Optional[str] = None
    structure: Optional[str] = None
    m_type: Optional[str] = None
    min_pos: int = 0
    max_pos: int = 999
    min_id: int = 0
    max_id: int = 999
    token: Optional[str] = None
    min_length: int = 0
    max_length: int = 999999
    strand_lengths: Optional[List[int]] = None
    has_children: Optional[bool] = None


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
            if msp.token is not None and m.token != msp.token:
                continue
            motif_length = m.end_pos - m.start_pos
            if motif_length < msp.min_length or motif_length > msp.max_length:
                continue
            if msp.strand_lengths is not None:
                lengths = [len(s) for s in m.strands]
                if lengths != msp.strand_lengths:
                    continue
            if msp.has_children is not None:
                has_children = len(m.children) > 0
                if has_children != msp.has_children:
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
            msp = MotifSearchParams(token=token)
        else:
            # Create a copy to avoid modifying the original
            msp = MotifSearchParams(
                sequence=msp.sequence,
                structure=msp.structure,
                m_type=msp.m_type,
                min_pos=msp.min_pos,
                max_pos=msp.max_pos,
                min_id=msp.min_id,
                max_id=msp.max_id,
                token=token,  # Override with provided token
                min_length=msp.min_length,
                max_length=msp.max_length,
                strand_lengths=msp.strand_lengths,
                has_children=msp.has_children,
            )
        return self.__get_motifs_by_params(msp)

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

    # Search methods #############################################################

    def find(self, sub: 'SecStruct', start: Optional[int] = None, end: Optional[int] = None) -> List[Tuple[int, int]]:
        """Find the position(s) of a substructure in this structure.

        Args:
            sub: The substructure to search for.
            start: Start position to search from (default: 0).
            end: End position to search to (default: end of sequence).

        Returns:
            List[Tuple[int, int]]: List of (start, end) tuples for matches.
                Each tuple represents the start (inclusive) and end (exclusive) positions.
        """
        if start is None:
            start = 0
        if end is None:
            end = len(self.__sequence)

        matches = []
        sub_seq = sub.__sequence
        sub_struct = sub.__structure

        # Handle multi-strand structures
        if "&" in sub_seq:
            # For multi-strand, need to match across strand boundaries
            # This is more complex - for now, just search within single strands
            search_seq = self.__sequence[start:end]
            search_struct = self.__structure[start:end]
            pos = search_seq.find(sub_seq)
            if pos != -1:
                # Verify structure matches at this position
                if search_struct[pos:pos + len(sub_seq)] == sub_struct:
                    matches.append((start + pos, start + pos + len(sub_seq)))
        else:
            # Single-strand search
            search_seq = self.__sequence[start:end]
            search_struct = self.__structure[start:end]
            pos = 0
            while True:
                pos = search_seq.find(sub_seq, pos)
                if pos == -1:
                    break
                # Verify structure matches at this position
                if search_struct[pos:pos + len(sub_seq)] == sub_struct:
                    matches.append((start + pos, start + pos + len(sub_seq)))
                pos += 1

        return matches

    def find_sequence(self, pattern: str, allow_wildcards: bool = True) -> List[Tuple[int, int]]:
        """Find positions matching a sequence pattern.

        Supports wildcards: 'N' matches any nucleotide, 'R' matches A/G, etc.

        Args:
            pattern: Sequence pattern to search for.
            allow_wildcards: If True, interpret wildcard characters (N, R, Y, etc.).

        Returns:
            List[Tuple[int, int]]: List of (start, end) tuples for matches.
        """
        matches = []
        if not allow_wildcards:
            # Simple string search
            pos = 0
            while True:
                pos = self.__sequence.find(pattern, pos)
                if pos == -1:
                    break
                matches.append((pos, pos + len(pattern)))
                pos += 1
        else:
            # Pattern matching with wildcards
            import re
            # Convert wildcards to regex
            regex_pattern = pattern.replace("N", "[AUCG]")
            regex_pattern = regex_pattern.replace("R", "[AG]")  # Purine
            regex_pattern = regex_pattern.replace("Y", "[UC]")  # Pyrimidine
            regex_pattern = regex_pattern.replace("M", "[AC]")  # Amino
            regex_pattern = regex_pattern.replace("K", "[UG]")  # Keto
            regex_pattern = regex_pattern.replace("S", "[GC]")  # Strong
            regex_pattern = regex_pattern.replace("W", "[AU]")  # Weak
            regex_pattern = regex_pattern.replace("B", "[UCG]")  # Not A
            regex_pattern = regex_pattern.replace("D", "[AUG]")  # Not C
            regex_pattern = regex_pattern.replace("H", "[AUC]")  # Not G
            regex_pattern = regex_pattern.replace("V", "[ACG]")  # Not U
            regex_pattern = regex_pattern.replace(".", "\\.")  # Escape dots

            for match in re.finditer(regex_pattern, self.__sequence):
                matches.append((match.start(), match.end()))

        return matches

    def find_structure(self, pattern: str) -> List[Tuple[int, int]]:
        """Find positions matching a structure pattern.

        Args:
            pattern: Structure pattern to search for (e.g., "(((...)))").

        Returns:
            List[Tuple[int, int]]: List of (start, end) tuples for matches.
        """
        matches = []
        pos = 0
        while True:
            pos = self.__structure.find(pattern, pos)
            if pos == -1:
                break
            matches.append((pos, pos + len(pattern)))
            pos += 1
        return matches

    # Connectivity and base pair methods #########################################

    @property
    def connectivity(self) -> List[int]:
        """Get connectivity list (pairmap).

        Returns:
            List[int]: Connectivity list where each index contains the paired
                position or -1 if unpaired.
        """
        from rna_secstruct.connectivity import connectivity_list
        return connectivity_list(self.__structure)

    def get_basepair(self, index: int) -> Optional[Tuple[int, int]]:
        """Get base pair for a given position.

        Args:
            index: The position index.

        Returns:
            Optional[Tuple[int, int]]: Tuple of (index, paired_index) if paired,
                None if unpaired.
        """
        conn = self.connectivity
        if conn[index] == -1:
            return None
        return (index, conn[index])

    def is_paired(self, index: int) -> bool:
        """Check if position is paired.

        Args:
            index: The position index.

        Returns:
            bool: True if the position is paired, False otherwise.
        """
        return self.connectivity[index] != -1

    # Statistics and analysis methods ############################################

    def get_num_basepairs(self) -> int:
        """Count number of base pairs.

        Returns:
            int: Number of base pairs in the structure.
        """
        conn = self.connectivity
        # Count pairs, but divide by 2 since each pair is counted twice
        return sum(1 for i, pair in enumerate(conn) if pair != -1 and i < pair)

    def get_num_unpaired(self) -> int:
        """Count number of unpaired nucleotides.

        Returns:
            int: Number of unpaired nucleotides.
        """
        conn = self.connectivity
        return sum(1 for pair in conn if pair == -1)

    def get_gc_content(self) -> float:
        """Calculate GC content.

        Returns:
            float: GC content as a fraction (0.0 to 1.0).
        """
        if len(self.__sequence) == 0:
            return 0.0
        gc_count = sum(1 for nuc in self.__sequence.upper() if nuc in 'GC')
        return gc_count / len(self.__sequence)

    def get_helix_lengths(self) -> List[int]:
        """Get lengths of all helices.

        Returns:
            List[int]: List of helix lengths.
        """
        helices = self.get_helices()
        return [len(h.sequence) for h in helices]

    # Validation utilities #######################################################

    def validate(self) -> None:
        """Validate structure is well-formed.

        Raises:
            ValueError: If structure is invalid with detailed explanation.
        """
        # Check length match
        if len(self.__sequence) != len(self.__structure):
            raise ValueError(
                f"Sequence and structure must have the same length. "
                f"Sequence length: {len(self.__sequence)}, structure length: {len(self.__structure)}."
            )

        # Check strand count match
        if self.__sequence.count("&") != self.__structure.count("&"):
            raise ValueError(
                f"Sequence and structure must have the same number of strands. "
                f"Sequence has {self.__sequence.count('&') + 1} strands, "
                f"structure has {self.__structure.count('&') + 1} strands."
            )

        # Validate structure format
        is_valid_dot_bracket_str(self.__structure)

        # Validate connectivity (this will raise if structure is malformed)
        try:
            _ = self.connectivity
        except Exception as e:
            raise ValueError(f"Invalid structure connectivity: {e}") from e

    def is_valid(self) -> bool:
        """Check if structure is valid (non-raising).

        Returns:
            bool: True if structure is valid, False otherwise.
        """
        try:
            self.validate()
            return True
        except (ValueError, TypeError):
            return False

    def normalize(self) -> 'SecStruct':
        """Return normalized version (uppercase, T->U conversion).

        Returns:
            SecStruct: New SecStruct instance with normalized sequence.
                Sequence is converted to uppercase and T is converted to U.
                Structure is unchanged.
        """
        normalized_seq = self.__sequence.upper().replace('T', 'U')
        return SecStruct(normalized_seq, self.__structure)

    # Comparison operations #######################################################

    def __eq__(self, other) -> bool:
        """Compare two structures for equality.

        Args:
            other: Another SecStruct to compare with.

        Returns:
            bool: True if sequences and structures are identical.
        """
        if not isinstance(other, SecStruct):
            return False
        return (
            self.__sequence == other.__sequence
            and self.__structure == other.__structure
        )

    def structural_similarity(self, other: 'SecStruct') -> float:
        """Calculate structural similarity score.

        Compares the structure strings and returns the fraction of positions
        that have the same structure character.

        Args:
            other: Another SecStruct to compare with.

        Returns:
            float: Similarity score between 0.0 and 1.0.
        """
        if len(self.__structure) != len(other.__structure):
            return 0.0
        if len(self.__structure) == 0:
            return 1.0

        matches = sum(
            1 for s1, s2 in zip(self.__structure, other.__structure) if s1 == s2
        )
        return matches / len(self.__structure)

    def sequence_identity(self, other: 'SecStruct') -> float:
        """Calculate sequence identity.

        Compares the sequences and returns the fraction of positions
        that have the same nucleotide.

        Args:
            other: Another SecStruct to compare with.

        Returns:
            float: Identity score between 0.0 and 1.0.
        """
        if len(self.__sequence) != len(other.__sequence):
            return 0.0
        if len(self.__sequence) == 0:
            return 1.0

        matches = sum(
            1 for s1, s2 in zip(self.__sequence, other.__sequence) if s1 == s2
        )
        return matches / len(self.__sequence)
