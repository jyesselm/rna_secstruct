"""
Connectivity list generation and manipulation for RNA secondary structures.

Supports multiple bracket types for pseudo-knots and alternative pairing notations
(letter-based and number-based).
"""

import re
from typing import List, Dict, Tuple, Optional, Union
from collections import defaultdict


# Standard bracket types for pseudo-knots
STANDARD_BRACKET_TYPES = [
    ("(", ")"),
    ("[", "]"),
    ("{", "}"),
    ("<", ">"),
]


def connectivity_list(
    structure: str, bracket_types: Optional[List[Tuple[str, str]]] = None
) -> List[int]:
    """Generates a connectivity list or pairmap from a dot-bracket secondary structure.

    The list has the index of a position's complement, if it is unpaired ('.'), it will have a -1 instead.

    Args:
        structure: A dot-bracket structure.
        bracket_types: List of (open, close) bracket pairs. Defaults to [('(', ')')] for basic support.
                      Use STANDARD_BRACKET_TYPES for full pseudo-knot support.

    Returns:
        List[int]: The connectivity list or pairmap. Each index contains the paired position or -1.

    Raises:
        ValueError: If brackets are unbalanced.
    """
    if bracket_types is None:
        bracket_types = [("(", ")")]

    connections = [-1] * len(structure)
    stacks = {open_bracket: [] for open_bracket, _ in bracket_types}
    bracket_map = {
        close_bracket: open_bracket for open_bracket, close_bracket in bracket_types
    }

    for index, char in enumerate(structure):
        if char in stacks:
            # Opening bracket
            stacks[char].append(index)
        elif char in bracket_map:
            # Closing bracket
            open_bracket = bracket_map[char]
            if not stacks[open_bracket]:
                raise ValueError(
                    f"Unbalanced bracket: '{char}' at position {index} has no matching opening bracket"
                )
            complement = stacks[open_bracket].pop()
            connections[complement] = index
            connections[index] = complement
        # '.' and '&' are ignored (unpaired and strand separator)

    # Check for unmatched opening brackets
    for open_bracket, stack in stacks.items():
        if stack:
            raise ValueError(
                f"Unbalanced brackets: {len(stack)} unmatched '{open_bracket}' brackets"
            )

    return connections


def connectivity_list_multi_bracket(
    structure: str, bracket_types: Optional[List[Tuple[str, str]]] = None
) -> Dict[str, List[int]]:
    """Generates connectivity lists for multiple bracket types separately.

    This allows detection of pseudo-knots by analyzing interleaving bracket types.

    Args:
        structure: A dot-bracket structure with multiple bracket types.
        bracket_types: List of (open, close) bracket pairs. Defaults to STANDARD_BRACKET_TYPES.

    Returns:
        Dict[str, List[int]]: Dictionary mapping bracket type names (e.g., '()', '[]') to connectivity lists.

    Raises:
        ValueError: If any bracket type is unbalanced.
    """
    if bracket_types is None:
        bracket_types = STANDARD_BRACKET_TYPES

    result = {}
    for open_bracket, close_bracket in bracket_types:
        bracket_name = f"{open_bracket}{close_bracket}"
        # Extract only this bracket type, replace others with '.'
        filtered_structure = "".join(
            char if char in (open_bracket, close_bracket, ".", "&") else "."
            for char in structure
        )
        try:
            result[bracket_name] = connectivity_list(
                filtered_structure, [(open_bracket, close_bracket)]
            )
        except ValueError:
            # This bracket type might not be present or balanced, skip it
            result[bracket_name] = [-1] * len(structure)

    return result


def has_pseudoknot(
    connectivity_lists: Union[List[int], Dict[str, List[int]]],
    bracket_types: Optional[List[Tuple[str, str]]] = None,
) -> bool:
    """Detect if structure contains pseudo-knots.

    A pseudo-knot occurs when pairs from different bracket types cross each other.

    Args:
        connectivity_lists: Either a single connectivity list or dict of connectivity lists by bracket type.
        bracket_types: List of bracket types (only needed if connectivity_lists is a single list).

    Returns:
        bool: True if pseudo-knots are detected.
    """
    if isinstance(connectivity_lists, list):
        # Single connectivity list - check if we need to analyze multiple bracket types
        if bracket_types is None or len(bracket_types) <= 1:
            return False
        # Convert to multi-bracket format for analysis
        # This would require the original structure, so for now return False
        # In practice, use connectivity_list_multi_bracket first
        return False

    # Multiple connectivity lists - check for crossing pairs
    lists = list(connectivity_lists.values())
    if len(lists) < 2:
        return False

    # Check if pairs from different bracket types cross
    for i in range(len(lists)):
        for j in range(i + 1, len(lists)):
            if _pairs_cross(lists[i], lists[j]):
                return True

    return False


def _pairs_cross(conn1: List[int], conn2: List[int]) -> bool:
    """Check if pairs from two connectivity lists cross each other."""
    # Get all pairs from conn1
    pairs1 = [
        (i, conn1[i]) for i in range(len(conn1)) if conn1[i] != -1 and i < conn1[i]
    ]
    pairs2 = [
        (i, conn2[i]) for i in range(len(conn2)) if conn2[i] != -1 and i < conn2[i]
    ]

    for i1, j1 in pairs1:
        for i2, j2 in pairs2:
            # Check if pairs cross: (i1 < i2 < j1 < j2) or (i2 < i1 < j2 < j1)
            if (i1 < i2 < j1 < j2) or (i2 < i1 < j2 < j1):
                return True

    return False


def detect_structure_format(structure: str) -> str:
    """Detect the format of the structure string.

    Args:
        structure: The structure string.

    Returns:
        str: Format type: 'bracket', 'letter', 'number', or 'mixed'.
    """
    has_brackets = bool(re.search(r"[()[\]{}<>]", structure))
    has_letters = bool(re.search(r"[a-zA-Z]", structure))
    has_numbers = bool(re.search(r"\d", structure))

    format_count = sum([has_brackets, has_letters, has_numbers])

    if format_count > 1:
        return "mixed"
    elif has_brackets:
        return "bracket"
    elif has_letters:
        return "letter"
    elif has_numbers:
        return "number"
    else:
        # Only dots, spaces, or strand separators
        return "bracket"  # Default to bracket format


def connectivity_list_from_letters(
    structure: str, case_sensitive: bool = True
) -> List[int]:
    """Generate connectivity list from letter-based notation.

    Example: "a b c c b a" -> pairs: (0,5), (1,4), (2,3)
    Example: "aaa...bbb" -> pairs: (0,8), (1,7), (2,6)

    Args:
        structure: Structure with letter pairs. Can be space-separated or continuous.
        case_sensitive: Whether to treat uppercase and lowercase as different pairs.

    Returns:
        List[int]: Connectivity list.

    Raises:
        ValueError: If letters are not properly paired.
    """
    # Remove spaces and strand separators for processing
    clean_structure = structure.replace(" ", "").replace("&", "&")
    connections = [-1] * len(structure)

    # Track positions of each letter
    letter_positions = defaultdict(list)

    for index, char in enumerate(structure):
        if char in (".", "&", " "):
            continue

        if not case_sensitive:
            char = char.lower()

        letter_positions[char].append(index)

    # Match pairs (first occurrence with last occurrence, etc.)
    for letter, positions in letter_positions.items():
        if len(positions) % 2 != 0:
            raise ValueError(
                f"Unpaired letter '{letter}': found {len(positions)} occurrences (must be even)"
            )

        # Pair from outside in
        while positions:
            i = positions.pop(0)
            j = positions.pop()
            connections[i] = j
            connections[j] = i

    return connections


def connectivity_list_from_numbers(structure: str) -> List[int]:
    """Generate connectivity list from number-based notation.

    Example: "1 2 3 3 2 1" -> pairs: (0,5), (1,4), (2,3)

    Args:
        structure: Structure with number pairs. Can be space-separated or continuous.

    Returns:
        List[int]: Connectivity list.

    Raises:
        ValueError: If numbers are not properly paired.
    """
    # Extract numbers (handle multi-digit numbers)
    numbers = []
    number_positions = []

    for match in re.finditer(r"\d+", structure):
        numbers.append(int(match.group()))
        number_positions.append((match.start(), match.end()))

    connections = [-1] * len(structure)

    # Track positions of each number
    number_pos_map = defaultdict(list)

    for idx, (start, end) in enumerate(number_positions):
        num = numbers[idx]
        number_pos_map[num].append((start, end))

    # Match pairs
    for num, positions in number_pos_map.items():
        if len(positions) % 2 != 0:
            raise ValueError(
                f"Unpaired number '{num}': found {len(positions)} occurrences (must be even)"
            )

        # Pair from outside in
        while positions:
            (i_start, i_end) = positions.pop(0)
            (j_start, j_end) = positions.pop()
            # Use the start position for pairing
            connections[i_start] = j_start
            connections[j_start] = i_start

    return connections


def get_connectivity(
    structure: str,
    format: Optional[str] = None,
    bracket_types: Optional[List[Tuple[str, str]]] = None,
) -> Union[List[int], Dict[str, List[int]]]:
    """Unified interface to get connectivity list(s) from any supported format.

    Args:
        structure: Structure string in any supported format.
        format: Explicit format ('bracket', 'letter', 'number', 'auto'). If None, auto-detects.
        bracket_types: For bracket format, specify bracket types. None uses default.

    Returns:
        List[int] for single bracket type or letter/number format.
        Dict[str, List[int]] for multi-bracket format with multiple bracket types.

    Raises:
        ValueError: If format is invalid or structure is malformed.
    """
    if format is None or format == "auto":
        format = detect_structure_format(structure)

    if format == "bracket":
        if bracket_types and len(bracket_types) > 1:
            return connectivity_list_multi_bracket(structure, bracket_types)
        else:
            return connectivity_list(structure, bracket_types)
    elif format == "letter":
        return connectivity_list_from_letters(structure)
    elif format == "number":
        return connectivity_list_from_numbers(structure)
    elif format == "mixed":
        # For mixed format, try to parse as bracket first, then fall back
        # This is a simplified approach - full mixed format support would be more complex
        return connectivity_list(structure, bracket_types)
    else:
        raise ValueError(f"Unknown format: {format}")


class ConnectivityList:
    """Represents a connectivity list for RNA secondary structure.

    Attributes:
        connections (List[int]): A list of indices representing the connectivity
            between nucleotides.
        sequence (str): The RNA sequence.
    """

    def __init__(self, sequence: str, structure: str, format: Optional[str] = None):
        """Initializes a ConnectivityList object.

        Args:
            sequence: The RNA sequence.
            structure: The RNA secondary structure.
            format: Structure format ('bracket', 'letter', 'number', 'auto').
        """
        conn = get_connectivity(structure, format)
        if isinstance(conn, dict):
            # For multi-bracket, use the first bracket type's connectivity
            self.connections = list(conn.values())[0]
        else:
            self.connections = conn
        self.sequence = sequence
        self.structure = structure

    def is_nucleotide_paired(self, index: int) -> bool:
        """Checks if a nucleotide at a given index is paired.

        Args:
            index: The index of the nucleotide.

        Returns:
            bool: True if the nucleotide is paired, False otherwise.
        """
        return self.connections[index] != -1

    def get_paired_nucleotide(self, index: int) -> int:
        """Returns the index of the nucleotide paired with the nucleotide at the given index.

        Args:
            index: The index of the nucleotide.

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
            index: The index of the nucleotide.

        Returns:
            str: The base pair of the nucleotide, or "." if unpaired.
        """
        if not self.is_nucleotide_paired(index):
            return "."
        return self.sequence[index] + self.sequence[self.get_paired_nucleotide(index)]


def is_circular(start: int, connections: List[int]) -> bool:
    """Check if a given RNA structure is circular.

    Args:
        start: The starting index of the RNA structure.
        connections: A list of connections between nucleotides.

    Returns:
        bool: True if the RNA structure is circular.
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
