"""
A simple parser of rna secondary structure inspired by `rna_library` code written
by Chris Jurich
"""

import re
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass

from rna_secstruct.motif import Motif
from rna_secstruct.logger import get_logger

log = get_logger("parser")


def is_valid_dot_bracket_str(structure: str) -> bool:
    """
    Checks if a structure is a valid dot-bracket structure containing only
    '(', '.' or ')' characters. Also checks if the structure is balanced.

    Args:
        structure (str): The dot bracket structure to be checked.

    Returns:
        bool: True if the structure is valid, False otherwise.

    Raises:
        ValueError: If the structure contains invalid characters or is unbalanced.
    """
    lparen_ct = 0
    for ch in structure:
        if ch == "(":
            lparen_ct += 1
        elif ch == ")":
            lparen_ct -= 1
        elif ch == "." or ch == "&":
            continue
        else:
            raise ValueError(
                f"{ch} is invalid in a dot-bracket structure. Only '(', '.' "
                f"and ')' are allowed"
            )

        if lparen_ct < 0:
            raise ValueError(f"{structure} is an unbalanced structure")

    if lparen_ct != 0:
        raise ValueError(f"{structure} is an unbalanced structure")

    for ii in range(3):
        invalid = "(" + "." * ii + ")"
        if structure.find(invalid) != -1:
            log.warning(f"{structure} has a hairpin that is too small")

    return True


# Import connectivity functions from connectivity module
from rna_secstruct.connectivity import (
    connectivity_list,
    ConnectivityList,
    is_circular,
    STANDARD_BRACKET_TYPES,
)

# Re-export for backward compatibility (deprecated - use connectivity module directly)
__all__ = ["Parser", "is_valid_dot_bracket_str", "connectivity_list", "ConnectivityList", "is_circular"]


def _detect_bracket_types(structure: str) -> List[Tuple[str, str]]:
    """Detect which bracket types are present in the structure.
    
    Args:
        structure: Structure string to analyze.
        
    Returns:
        List of (open, close) bracket pairs found in structure.
        Defaults to [('(', ')')] if no brackets detected.
    """
    detected = []
    for open_b, close_b in STANDARD_BRACKET_TYPES:
        if open_b in structure or close_b in structure:
            detected.append((open_b, close_b))
    
    # Default to standard parentheses if no brackets detected
    return detected if detected else [("(", ")")]


class Parser:
    """A class to parse secondary structure into motifs."""

    def __init__(self):
        self.motif_id = 0

    def parse(self, sequence: str, structure: str) -> None:
        """
        Parse the given sequence and structure into motifs.

        Args:
            sequence: A sequence of nucleotides.
            structure: A dot bracket structure.

        Returns:
            None
        """
        self.motif_id = 0
        self.__check_to_see_if_inputs_valid(sequence, structure)
        # Detect bracket types and pass to connectivity_list
        # For single bracket type (standard parentheses), pass None for backward compatibility
        # For multiple bracket types, pass the list
        bracket_types = _detect_bracket_types(structure)
        if len(bracket_types) == 1 and bracket_types[0] == ("(", ")"):
            # Default behavior for backward compatibility
            connections = connectivity_list(structure, bracket_types=None)
        else:
            # Multiple bracket types or non-standard brackets
            connections = connectivity_list(structure, bracket_types=bracket_types)
        return self.__get_motifs(sequence, structure, connections, 0)

    def __check_to_see_if_inputs_valid(self, sequence: str, structure: str) -> None:
        """
        Check if the inputs are valid.

        Args:
            sequence: A sequence of nucleotides.
            structure: A dot bracket structure.

        Returns:
            None

        Raises:
            ValueError: If the sequence is empty or if the sequence and structure
                have different lengths.
        """
        if len(sequence) == 0:
            raise ValueError("Sequence is empty")
        # enforce upper case sequence and is RNA
        sequence = sequence.upper().replace("T", "U")
        if len(sequence) != len(structure):
            raise ValueError(
                f"sequence and structure are not the same length: {sequence} {structure}"
            )
        if not re.match(r"^[ACGUTN&]+$", sequence):
            log.warning(f"sequence contains invalid characters: {sequence}")
        
        # Allow standard bracket types: (), [], {}, <>
        allowed_chars = set(".& ")
        for open_b, close_b in STANDARD_BRACKET_TYPES:
            allowed_chars.add(open_b)
            allowed_chars.add(close_b)
        
        # Check if structure contains only allowed characters
        invalid_chars = [c for c in structure if c not in allowed_chars]
        if invalid_chars:
            raise ValueError(
                f"structure contains invalid characters: {structure}. "
                f"Invalid characters found: {set(invalid_chars)}. "
                f"Allowed characters: {sorted(allowed_chars)}"
            )
        
        # Note: Bracket balance validation is done in parse() method via connectivity_list()
        # This avoids calling connectivity_list() twice (once here, once in parse())

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
        motifs = []
        if start >= len(connections):
            return None
        if connections[start] == -1:
            return self.__get_single_strand(sequence, structure, connections, start)
        return self.__get_helix(sequence, structure, connections, start)

    def __get_single_strand(
        self, sequence: str, structure: str, connections: List[int], start: int
    ) -> Motif:
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
                self.__get_motifs(
                    sequence, structure, connections, start + single_strand_count
                )
            )
        return sstrand

    def __get_helix(
        self, sequence: str, structure: str, connections: List[int], start: int
    ) -> Motif:
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
        helix_len = self.__get_helix_length(connections, start)
        if helix_len == 0:
            # This should not happen if structure is valid, but handle gracefully
            raise ValueError(f"No valid helix found starting at position {start}. Structure may be invalid.")
        lhs, rhs = [], []
        for index in range(start, start + helix_len):
            lhs.append(index)
            rhs.append(connections[index])
        rhs.reverse()
        seq1, ss1 = self.__get_seq_and_ss_from_strand(sequence, structure, lhs)
        seq2, ss2 = self.__get_seq_and_ss_from_strand(sequence, structure, rhs)
        helix = Motif(
            "HELIX", [lhs, rhs], f"{seq1}&{seq2}", f"{ss1}&{ss2}", self.motif_id
        )
        self.motif_id += 1
        # Check if there's a junction/hairpin inside the helix
        # The junction/hairpin starts at the last position of the helix (the inner base pair)
        if helix_len > 0 and start + helix_len - 1 < len(connections):
            if connections[start + helix_len - 1] > start:
                junction_start = start + helix_len - 1
                # Add bounds check to ensure junction_start is valid
                if junction_start < len(connections):
                    helix.add_child(
                        self.__get_junction_or_hairpin(
                            sequence, structure, connections, junction_start
                        )
                    )

        if len(rhs) > 0 and not is_circular(rhs[-1], connections):
            next_start = rhs[-1] + 1
            if next_start < len(connections) and next_start > start:
                motif = self.__get_motifs(sequence, structure, connections, next_start)
                if motif is not None:
                    helix.add_child(motif)
        return helix

    def __get_junction_or_hairpin(
        self, sequence: str, structure: str, connections: List[int], start: int
    ) -> Motif:
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
        if start < 0 or start >= len(structure):
            raise ValueError(f"start position {start} is out of bounds")
        strands = []
        pos = start
        # pos should be the first opening pair of a junction or hairpin
        # Support any opening bracket from STANDARD_BRACKET_TYPES
        opening_brackets = {open_b for open_b, _ in STANDARD_BRACKET_TYPES}
        if structure[pos] not in opening_brackets:
            raise ValueError(
                f"expected opening bracket at position {pos}, got '{structure[pos]}'. "
                f"Valid opening brackets: {sorted(opening_brackets)}"
            )
        iterations = 0
        max_iterations = len(connections)  # Prevent infinite loops
        while True:
            if iterations >= max_iterations:
                raise ValueError(f"infinite loop detected in junction/hairpin parsing at position {start}")
            iterations += 1
            if pos < 0 or pos >= len(connections):
                raise ValueError(f"position {pos} is out of bounds")
            next_strand = [pos]
            pos += 1
            if pos >= len(connections):
                raise ValueError(f"position {pos} is out of bounds")
            while pos < len(connections) and connections[pos] == -1:
                next_strand.append(pos)
                pos += 1
            if pos >= len(connections):
                raise ValueError(f"position {pos} is out of bounds")
            next_strand.append(pos)
            strands.append(next_strand)
            if connections[pos] < 0 or connections[pos] >= len(connections):
                raise ValueError(f"connection at position {pos} is invalid: {connections[pos]}")
            pos = connections[pos]
            # made a complete circle
            if pos == start:
                break
        # is a junction
        self.motif_id += 1
        if len(strands) > 1:
            seq_and_ss = [
                self.__get_seq_and_ss_from_strand(sequence, structure, strand)
                for strand in strands
            ]
            seq = "&".join([seq for seq, ss in seq_and_ss])
            ss = "&".join([ss for seq, ss in seq_and_ss])
            m = Motif("JUNCTION", strands, seq, ss, self.motif_id - 1)
            for strand in strands[:-1]:
                next_pos = strand[-1] + 1
                if next_pos < len(connections) and next_pos > start:
                    child = self.__get_motifs(sequence, structure, connections, next_pos)
                    if child is not None:
                        m.add_child(child)
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
        if start >= len(connections) or start < 0:
            return 0
        complement = connections[start]
        if complement < 0 or complement >= len(connections):
            return 0
        length = 0
        max_len = len(connections)
        while (
            start + length < max_len
            and complement - length >= 0
            and connections[start + length] == complement - length
            and connections[complement - length] == start + length
        ):
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
        seq = "".join([sequence[i] for i in strand])
        ss = "".join([structure[i] for i in strand])
        return seq, ss
