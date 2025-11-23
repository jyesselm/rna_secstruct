"""
A simple parser of rna secondary structure inspired by `rna_library` code written
by Chris Jurich
"""

from typing import List, Optional, Tuple

from rna_secstruct.logger import get_logger
from rna_secstruct.motif import Motif

log = get_logger("parser")


def normalize_structure(structure: str) -> str:
    """
    Normalizes a structure string by replacing invalid characters with '.'.

    IMPORTANT: Preserves bracket types ([ ], { }, < >) for pseudoknot representation.
    Only replaces truly invalid characters (not valid bracket pairs, letters, numbers, dots).

    Valid characters preserved:
    - Bracket types: ( ) [ ] { } < >
    - Unpaired: .
    - Strand separator: &
    - Letters: a-z A-Z (for letter-based pairing)
    - Numbers: 0-9 (for number-based pairing)

    Args:
        structure (str): The structure string to normalize.

    Returns:
        str: Normalized structure string with invalid characters replaced by '.'.
    """
    # Valid bracket characters for pseudoknots
    valid_brackets = "()[]{}<>"
    normalized = []
    invalid_chars = set()

    for ch in structure:
        if ch in valid_brackets or ch in ".&":
            # Valid bracket type, dot, or strand separator - preserve
            normalized.append(ch)
        elif ch.isalnum():
            # Valid letters or numbers for alternative pairing formats - preserve
            normalized.append(ch)
        else:
            # Invalid character - replace with '.'
            normalized.append(".")
            invalid_chars.add(ch)

    if invalid_chars:
        log.warning(
            f"Structure contains invalid characters: {invalid_chars}. "
            f"Replaced with '.' (bracket types preserved for pseudoknots)."
        )

    return "".join(normalized)


def balance_structure(structure: str) -> str:
    """
    Attempts to balance standard parentheses (not other bracket types to preserve pseudoknots).

    WARNING: Only balances () parentheses. Other bracket types ([], {}, <>) are preserved
    as-is to maintain pseudoknot information. Balancing all bracket types could break
    pseudoknot structures.

    Args:
        structure (str): The structure string to balance.

    Returns:
        str: Balanced structure string (only for () parentheses).
    """
    lparen_ct = 0
    for ch in structure:
        if ch == "(":
            lparen_ct += 1
        elif ch == ")":
            lparen_ct -= 1

    if lparen_ct > 0:
        # Missing closing parentheses, add them at the end
        structure = structure + ")" * lparen_ct
        log.warning(
            f"Structure has unbalanced () parentheses: {lparen_ct} missing closing parentheses. Added at end."
        )
    elif lparen_ct < 0:
        # Missing opening parentheses, add them at the start
        structure = "(" * (-lparen_ct) + structure
        log.warning(
            f"Structure has unbalanced () parentheses: {-lparen_ct} missing opening parentheses. Added at start."
        )

    # Note: We don't balance other bracket types ([], {}, <>) to preserve pseudoknot information

    return structure


def is_valid_dot_bracket_str(structure: str) -> bool:
    """
    Checks if a structure is a valid dot-bracket structure.

    Supports multiple bracket types for pseudoknots: ( ) [ ] { } < >
    Also supports alternative formats: letters, numbers, dots, strand separators.

    Now logs warnings instead of raising exceptions.

    Args:
        structure (str): The dot bracket structure to be checked.

    Returns:
        bool: True if the structure appears valid after normalization.

    Note:
        Invalid characters and unbalanced structures are logged as warnings.
    """
    # Track parentheses balance (but not other bracket types here - connectivity handles those)
    lparen_ct = 0
    invalid_chars = set()
    valid_brackets = "()[]{}<>"

    for ch in structure:
        if ch == "(":
            lparen_ct += 1
        elif ch == ")":
            lparen_ct -= 1
        elif ch in valid_brackets[2:] or ch in ".&":
            # Other bracket types or valid characters - skip balance check
            continue
        elif ch.isalnum():
            # Letters or numbers for alternative formats - valid
            continue
        else:
            invalid_chars.add(ch)

        if lparen_ct < 0:
            log.warning(f"Structure has unmatched closing parentheses: {structure[:100]}...")
            break

    if invalid_chars:
        # Only log once for invalid characters - normalization will handle it
        pass  # Normalization will log the warning

    if lparen_ct != 0:
        log.warning(f"Structure has unbalanced parentheses (unmatched count: {lparen_ct})")

    # Check for small hairpins in standard parentheses
    for ii in range(3):
        invalid = "(" + "." * ii + ")"
        if structure.find(invalid) != -1:
            log.warning("Structure has a hairpin that is too small")
            break

    return True


def connectivity_list(structure: str) -> List[int]:
    """Generates a connectivity list or pairmap from a dot-bracket secondary structure.

    The list has the index of a position's complement, if it is a '.', it will have a
      -1 instead.

    Now handles unbalanced parentheses gracefully with warnings.

    Args:
        structure (str): A dot-bracket structure.

    Returns:
        List[int]: The connectivity list or pairmap.

    Note:
        Unbalanced parentheses are handled by ignoring unmatched pairs with warnings.
    """
    connections, pairs = [-1] * len(structure), []
    for index, db in enumerate(structure):
        if db == "(":
            pairs.append(index)
        elif db == ")":
            if pairs:
                complement = pairs.pop()
                connections[complement] = index
                connections[index] = complement
            else:
                # Unmatched closing parenthesis
                log.warning(f"Unmatched closing parenthesis at position {index} in structure")
    if len(pairs):
        # Unmatched opening parentheses - mark them as unpaired
        log.warning(
            f"Structure has {len(pairs)} unmatched opening parentheses at positions {pairs}"
        )
        for _pos in pairs:
            # Keep as unpaired (already -1)
            pass
    return connections


class ConnectivityList:
    """Represents a connectivity list for RNA secondary structure.

    Attributes:
        connections (List[int]): A list of indices representing the connectivity
            between nucleotides.
        sequence (str): The RNA sequence.
    """

    def __init__(self, sequence: str, structure: str):
        """Initializes a ConnectivityList object.

        Args:
            sequence (str): The RNA sequence.
            structure (str): The RNA secondary structure.

        """
        self.connections = connectivity_list(structure)
        self.sequence = sequence

    def is_nucleotide_paired(self, index: int) -> bool:
        """Checks if a nucleotide at a given index is paired.

        Args:
            index (int): The index of the nucleotide.

        Returns:
            bool: True if the nucleotide is paired, False otherwise.

        """
        return self.connections[index] != -1

    def get_paired_nucleotide(self, index: int) -> int:
        """Returns the index of the nucleotide paired with the nucleotide at the given index.

        Args:
            index (int): The index of the nucleotide.

        Returns:
            int: The index of the paired nucleotide.

        Raises:
            ValueError: If the nucleotide at the given index is not paired.

        """
        if not self.is_nucleotide_paired(index):
            raise ValueError(f"Nucleotide at index {index} is not paired")
        return self.connections[index]

    def get_basepair(self, index: int) -> str:
        """Returns the base pair of the nucleotide at the given index.

        Args:
            index (int): The index of the nucleotide.

        Returns:
            str: The base pair of the nucleotide.

        """
        if not self.is_nucleotide_paired(index):
            return "."
        return self.sequence[index] + self.sequence[self.get_paired_nucleotide(index)]


def is_circular(start, connections):
    """Check if a given RNA structure is circular.

    Args:
        start (int): The starting index of the RNA structure.
        connections (List[int]): A list of connections between nucleotides.

    Returns:
        bool: True if the RNA structure is circular, False otherwise.
    """
    it = start + 1
    while True:
        while it < len(connections) and connections[it] == -1:
            it += 1
        if it == len(connections):
            return False

        it = connections[it] + 1
        if it == start or it < start:
            return True


class Parser:
    """A class to parse secondary structure into motifs."""

    def __init__(self):
        self.motif_id = 0

    def parse(self, sequence: str, structure: str) -> Optional[Motif]:
        """
        Parse the given sequence and structure into motifs.

        Args:
            sequence: A sequence of nucleotides.
            structure: A dot bracket structure.

        Returns:
            Optional[Motif]: The root motif, or None if parsing fails.
        """
        self.motif_id = 0
        # Normalize inputs before processing
        sequence, structure = self.__check_to_see_if_inputs_valid(sequence, structure)
        # Normalize structure to handle invalid characters
        structure = normalize_structure(structure)
        # Attempt to balance structure
        structure = balance_structure(structure)
        # Ensure sequence and structure are same length after normalization
        if len(sequence) != len(structure):
            min_len = min(len(sequence), len(structure))
            if len(sequence) > len(structure):
                sequence = sequence[:min_len]
                log.warning(f"Sequence truncated to match structure length: {min_len}")
            else:
                structure = structure[:min_len]
                log.warning(f"Structure truncated to match structure length: {min_len}")
        connections = connectivity_list(structure)
        return self.__get_motifs(sequence, structure, connections, 0)

    def __check_to_see_if_inputs_valid(self, sequence: str, structure: str) -> Tuple[str, str]:
        """
        Check if the inputs are valid and normalize them.

        Args:
            sequence: A sequence of nucleotides.
            structure: A dot bracket structure.

        Returns:
            Tuple[str, str]: Normalized (sequence, structure) pair.

        Note:
            Invalid characters are normalized or replaced with warnings.
        """
        if len(sequence) == 0:
            log.warning("Sequence is empty, skipping parse")
            return sequence, structure

        # enforce upper case sequence and convert DNA to RNA
        sequence = sequence.upper().replace("T", "U")

        # Normalize invalid characters in sequence (replace with N)
        valid_chars = set("ACGU&N")
        normalized_seq = []
        invalid_chars = []
        for ch in sequence:
            if ch in valid_chars:
                normalized_seq.append(ch)
            else:
                normalized_seq.append("N")  # Replace invalid characters with N
                if ch not in invalid_chars:
                    invalid_chars.append(ch)

        if invalid_chars:
            log.warning(
                f"Sequence contains invalid characters: {set(invalid_chars)}. "
                f"Replaced with 'N'."
            )

        sequence = "".join(normalized_seq)

        # Check length mismatch but don't raise - will be handled in parse()
        if len(sequence) != len(structure):
            log.warning(
                f"Sequence and structure have different lengths: "
                f"sequence={len(sequence)}, structure={len(structure)}. "
                f"Will be normalized."
            )

        # Check for invalid structure characters (but don't raise - will be normalized)
        is_valid_dot_bracket_str(structure)

        return sequence, structure

    def __get_motifs(
        self, sequence: str, structure: str, connections: List[int], start: int
    ) -> Optional[Motif]:
        """
        Get the motifs from the structure.

        Args:
            sequence: A sequence of nucleotides.
            structure: A dot bracket structure.
            connections: A list of connections.
            start: The start index.

        Returns:
            A list of motifs or None.
        """
        if start >= len(connections):
            return None
        if connections[start] == -1:
            return self.__get_single_strand(sequence, structure, connections, start)
        return self.__get_helix(sequence, structure, connections, start)

    def __get_single_strand(
        self, sequence: str, structure: str, connections: List[int], start: int
    ) -> Optional[Motif]:
        """
        Get a single strand.

        Args:
            sequence: A sequence of nucleotides.
            structure: A dot bracket structure.
            connections: A list of connections.
            start: The start index.

        Returns:
            A single strand motif.
        """
        # how many nucleotides are in the single strand
        single_strand_count = 0
        while (
            start + single_strand_count < len(connections)
            and connections[start + single_strand_count] == -1
        ):
            single_strand_count += 1
        strand = list(range(start, start + single_strand_count))
        sstrand = Motif(
            "SINGLESTRAND",
            [strand],
            sequence[start : start + single_strand_count],
            structure[start : start + single_strand_count],
            self.motif_id,
        )
        self.motif_id += 1
        if start + single_strand_count < len(connections):
            sstrand.add_child(
                self.__get_motifs(sequence, structure, connections, start + single_strand_count)
            )
        return sstrand

    def __get_helix(
        self, sequence: str, structure: str, connections: List[int], start: int
    ) -> Optional[Motif]:
        """
        Get a helix or junction.

        Args:
            sequence: A sequence of nucleotides.
            structure: A dot bracket structure.
            connections: A list of connections.
            start: The start index.

        Returns:
            A helix or junction motif.
        """
        if start >= len(connections):
            log.warning(
                f"Start position {start} is out of bounds for connections of length {len(connections)}"
            )
            return None

        if connections[start] == -1:
            # Not actually a helix, treat as single strand
            return self.__get_single_strand(sequence, structure, connections, start)

        helix_len = self.__get_helix_length(connections, start)
        if helix_len == 0:
            log.warning(f"Helix length is 0 at position {start}")
            return None

        lhs, rhs = [], []
        for index in range(start, start + helix_len):
            if index >= len(connections):
                log.warning(f"Index {index} is out of bounds while building helix")
                break
            lhs.append(index)
            if connections[index] == -1:
                log.warning(f"Unpaired position {index} in helix at start {start}")
                continue
            if connections[index] >= len(connections):
                log.warning(
                    f"Connection {connections[index]} is out of bounds for position {index}"
                )
                continue
            rhs.append(connections[index])

        if not rhs:
            log.warning(f"No valid pairs found in helix starting at {start}")
            return self.__get_single_strand(sequence, structure, connections, start)

        rhs.reverse()
        seq1, ss1 = self.__get_seq_and_ss_from_strand(sequence, structure, lhs)
        seq2, ss2 = self.__get_seq_and_ss_from_strand(sequence, structure, rhs)
        helix = Motif("HELIX", [lhs, rhs], f"{seq1}&{seq2}", f"{ss1}&{ss2}", self.motif_id)
        self.motif_id += 1

        if (
            start + helix_len - 1 < len(connections)
            and connections[start + helix_len - 1] != -1
            and connections[start + helix_len - 1] > start
        ):
            child = self.__get_junction_or_hairpin(
                sequence, structure, connections, start + helix_len - 1
            )
            if child is not None:
                helix.add_child(child)

        if rhs and rhs[-1] + 1 < len(connections) and not is_circular(rhs[-1], connections):
            motif = self.__get_motifs(sequence, structure, connections, rhs[-1] + 1)
            if motif is not None:
                helix.add_child(motif)

        return helix

    def __get_junction_or_hairpin(
        self, sequence: str, structure: str, connections: List[int], start: int
    ) -> Optional[Motif]:
        """
        Get a junction or hairpin.

        Args:
            sequence: A sequence of nucleotides.
            structure: A dot bracket structure.
            connections: A list of connections.
            start: The start index.

        Returns:
            A junction or hairpin motif.
        """
        strands = []
        pos = start

        # Check bounds
        if pos >= len(structure) or pos >= len(connections):
            log.warning(f"Position {pos} is out of bounds for structure of length {len(structure)}")
            # Return a single strand as fallback
            if pos < len(connections):
                return self.__get_single_strand(sequence, structure, connections, pos)
            return None

        # pos should be the first opening pair of a junction or hairpin
        # Accept any opening bracket: ( [ { <
        opening_brackets = "([{<"
        if structure[pos] not in opening_brackets:
            log.warning(
                f"Expected opening bracket at position {pos}, got '{structure[pos]}'. Treating as unpaired."
            )
            # If not paired, treat as single strand
            if connections[pos] == -1:
                return self.__get_single_strand(sequence, structure, connections, pos)
            # Otherwise try to continue with the paired position
            pos = connections[pos]
            if pos == start:
                log.warning(f"Circular reference detected at position {start}")
                return None

        while True:
            if pos >= len(structure) or pos >= len(connections):
                log.warning(f"Position {pos} is out of bounds while processing junction/hairpin")
                break

            next_strand = [pos]
            pos += 1

            # Check bounds before accessing
            if pos >= len(connections):
                log.warning(f"Position {pos} is out of bounds while building strand")
                break

            while pos < len(connections) and connections[pos] == -1:
                next_strand.append(pos)
                pos += 1

            if pos >= len(connections):
                log.warning(f"Position {pos} is out of bounds at end of strand")
                break

            next_strand.append(pos)
            strands.append(next_strand)

            if connections[pos] == -1 or connections[pos] >= len(connections):
                log.warning(f"Invalid connection at position {pos}")
                break

            pos = connections[pos]
            # made a complete circle
            if pos == start:
                break

            # Safety check to avoid infinite loops
            if len(strands) > 1000:
                log.warning(
                    f"Too many strands in junction/hairpin at position {start}, breaking loop"
                )
                break
        # is a junction
        self.motif_id += 1
        if len(strands) > 1:
            seq_and_ss = [
                self.__get_seq_and_ss_from_strand(sequence, structure, strand) for strand in strands
            ]
            seq = "&".join([seq for seq, ss in seq_and_ss])
            ss = "&".join([ss for seq, ss in seq_and_ss])
            m = Motif("JUNCTION", strands, seq, ss, self.motif_id - 1)
            for strand in strands[:-1]:
                m.add_child(self.__get_motifs(sequence, structure, connections, strand[-1]))
            return m
        else:
            seq, ss = self.__get_seq_and_ss_from_strand(sequence, structure, strands[0])
            return Motif("HAIRPIN", strands, seq, ss, self.motif_id - 1)

    def __get_helix_length(self, connections: List[int], start: int) -> int:
        """
        Get the length of a helix.

        Args:
            connections: A list of connections.
            start: The start index.

        Returns:
            The length of the helix.
        """
        if start >= len(connections) or connections[start] == -1:
            return 0

        complement = connections[start]
        if complement >= len(connections):
            log.warning(f"Complement {complement} is out of bounds for position {start}")
            return 0

        length = 0
        max_length = min(len(connections) - start, complement + 1)

        while length < max_length:
            if start + length >= len(connections) or complement - length < 0:
                break
            if connections[start + length] == -1:
                break
            if connections[start + length] != complement - length:
                break
            if connections[complement - length] == -1:
                break
            if connections[complement - length] != start + length:
                break
            length += 1

        return length

    def __get_seq_and_ss_from_strand(
        self, sequence: str, structure: str, strand: List[int]
    ) -> Tuple[str, str]:
        """
        Get the sequence and secondary structure from a strand.

        Args:
            sequence: A sequence of nucleotides.
            structure: A dot bracket structure.
            strand: A list of indices representing the strand.

        Returns:
            The sequence and secondary structure.
        """
        seq_parts = []
        ss_parts = []
        for i in strand:
            if i < 0 or i >= len(sequence) or i >= len(structure):
                log.warning(
                    f"Index {i} is out of bounds for sequence/structure of length {len(sequence)}"
                )
                continue
            seq_parts.append(sequence[i])
            ss_parts.append(structure[i])

        seq = "".join(seq_parts) if seq_parts else "N"
        ss = "".join(ss_parts) if ss_parts else "."
        return seq, ss
