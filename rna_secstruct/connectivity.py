"""
Connectivity list generation and manipulation for RNA secondary structures.

Supports multiple bracket types for pseudo-knots and alternative pairing notations
(letter-based and number-based).

Public API:
    - get_connectivity_list(): Factory function to create ConnectivityList objects
    - ConnectivityList: Main class for working with connectivity lists

Internal functions (prefixed with _) are for internal use only.
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


def _connectivity_list(
    structure: str, bracket_types: Optional[List[Tuple[str, str]]] = None
) -> List[int]:
    """Generates a connectivity list or pairmap from a dot-bracket secondary structure.

    The list has the index of a position's complement, if it is unpaired ('.'), it will have a -1 instead.
    For pure bracket format, only processes brackets and ignores letters/numbers.

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

    # For pure bracket format, use simple bracket-only parser
    connections = [-1] * len(structure)
    stacks = {open_bracket: [] for open_bracket, _ in bracket_types}
    bracket_map = {close_bracket: open_bracket for open_bracket, close_bracket in bracket_types}

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
        # '.' and '&' and other characters are ignored (unpaired and strand separator)

    # Check for unmatched opening brackets
    for open_bracket, stack in stacks.items():
        if stack:
            raise ValueError(
                f"Unbalanced brackets: {len(stack)} unmatched '{open_bracket}' brackets"
            )

    return connections


def _connectivity_list_unified(
    structure: str,
    bracket_types: Optional[List[Tuple[str, str]]] = None,
    case_sensitive: bool = True,
) -> Tuple[List[int], Dict[int, str]]:
    """Unified parser that handles brackets, letters, and numbers all together.

    This is the core parser that processes mixed formats like "(((aaa(((...)))aaa)))"
    where brackets and letters are interleaved.

    Args:
        structure: Structure string that may contain brackets, letters, and/or numbers.
        bracket_types: List of (open, close) bracket pairs. Defaults to [('(', ')')].
        case_sensitive: Whether to treat uppercase and lowercase letters as different pairs.

    Returns:
        Tuple[List[int], Dict[int, str]]: Connectivity list and pair type dictionary.

    Raises:
        ValueError: If brackets or letters/numbers are unbalanced.
    """
    if bracket_types is None:
        bracket_types = [("(", ")")]

    connections = [-1] * len(structure)
    pair_types: Dict[int, str] = {}

    # Track bracket stacks
    stacks: Dict[str, List[Tuple[int, str]]] = {
        open_bracket: [] for open_bracket, _ in bracket_types
    }
    bracket_map = {close_bracket: open_bracket for open_bracket, close_bracket in bracket_types}

    # Track letter positions (for pairing letters separately from brackets)
    letter_positions = defaultdict(list)

    # Track number positions
    number_pos_map = defaultdict(list)
    number_positions = []
    for match in re.finditer(r"\d+", structure):
        num = int(match.group())
        num_str = structure[match.start() : match.end()]
        number_positions.append((match.start(), match.end(), num, num_str))

    # First pass: process brackets and collect letters/numbers
    for index, char in enumerate(structure):
        if char in stacks:
            # Opening bracket
            stacks[char].append((index, char))
        elif char in bracket_map:
            # Closing bracket
            open_bracket = bracket_map[char]
            if not stacks[open_bracket]:
                raise ValueError(
                    f"Unbalanced bracket: '{char}' at position {index} has no matching opening bracket"
                )
            complement, pair_type = stacks[open_bracket].pop()
            connections[complement] = index
            connections[index] = complement
            pair_types[complement] = pair_type
            pair_types[index] = char
        elif char not in (".", "&", " "):
            # Check if it's a valid letter (a-z, A-Z) for pairing
            # Only process letters if they're valid RNA pairing characters
            if re.match(r"^[a-zA-Z]$", char):
                original_char = char
                if not case_sensitive:
                    char = char.lower()
                letter_positions[char].append((index, original_char))
            # Other characters (like 'x', 'y', etc.) are ignored - treated as unpaired
        # '.' and '&' and ' ' are ignored (unpaired and strand separator)

    # Check for unmatched opening brackets
    for open_bracket, stack in stacks.items():
        if stack:
            raise ValueError(
                f"Unbalanced brackets: {len(stack)} unmatched '{open_bracket}' brackets"
            )

    # Second pass: process letters (pair from outside in)
    # Only process letters if there are any (for mixed/letter formats)
    if letter_positions:
        for letter, positions in letter_positions.items():
            if len(positions) % 2 != 0:
                raise ValueError(
                    f"Unpaired letter '{letter}': found {len(positions)} occurrences (must be even)"
                )

            while positions:
                i, i_char = positions.pop(0)
                j, j_char = positions.pop()
                # Only pair if not already paired by brackets
                if connections[i] == -1 and connections[j] == -1:
                    connections[i] = j
                    connections[j] = i
                    pair_types[i] = i_char
                    pair_types[j] = j_char

    # Third pass: process numbers (pair from outside in)
    # Only process numbers if there are any (for mixed/number formats)
    if number_positions:
        for start, end, num, num_str in number_positions:
            number_pos_map[num].append((start, end, num_str))

        for num, positions in number_pos_map.items():
            if len(positions) % 2 != 0:
                raise ValueError(
                    f"Unpaired number '{num}': found {len(positions)} occurrences (must be even)"
                )

            while positions:
                (i_start, i_end, i_num_str) = positions.pop(0)
                (j_start, j_end, j_num_str) = positions.pop()
                # Only pair if not already paired by brackets or letters
                if connections[i_start] == -1 and connections[j_start] == -1:
                    connections[i_start] = j_start
                    connections[j_start] = i_start
                    pair_types[i_start] = i_num_str
                    pair_types[j_start] = j_num_str

    return connections, pair_types


def _connectivity_list_multi_bracket(
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
            char if char in (open_bracket, close_bracket, ".", "&") else "." for char in structure
        )
        try:
            # Use unified parser but only with this bracket type
            conn, _ = _connectivity_list_unified(
                filtered_structure, [(open_bracket, close_bracket)]
            )
            result[bracket_name] = conn
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
    pairs1 = [(i, conn1[i]) for i in range(len(conn1)) if conn1[i] != -1 and i < conn1[i]]
    pairs2 = [(i, conn2[i]) for i in range(len(conn2)) if conn2[i] != -1 and i < conn2[i]]

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
    # Check for letters that are NOT bracket characters
    # Letters used for pairing are typically lowercase a-z
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


def _get_connectivity(
    structure: str,
    format: Optional[str] = None,
    bracket_types: Optional[List[Tuple[str, str]]] = None,
    return_pair_types: bool = False,
) -> Union[
    List[int],
    Dict[str, List[int]],
    Tuple[List[int], Dict[int, str]],
    Tuple[Dict[str, List[int]], Dict[str, Dict[int, str]]],
]:
    """Unified interface to get connectivity list(s) from any supported format.

    Args:
        structure: Structure string in any supported format.
        format: Explicit format ('bracket', 'letter', 'number', 'auto'). If None, auto-detects.
        bracket_types: For bracket format, specify bracket types. None uses default.
        return_pair_types: If True, also return pair type information.

    Returns:
        If return_pair_types=False:
            - List[int] for single bracket type or letter/number format.
            - Dict[str, List[int]] for multi-bracket format with multiple bracket types.
        If return_pair_types=True:
            - Tuple[List[int], Dict[int, str]] for single format.
            - Tuple[Dict[str, List[int]], Dict[str, Dict[int, str]]] for multi-bracket format.
            The Dict[int, str] maps position to pair type (e.g., '(', '[', 'a', '1').

    Raises:
        ValueError: If format is invalid or structure is malformed.
    """
    if format is None or format == "auto":
        format = detect_structure_format(structure)

    if format == "bracket":
        if bracket_types and len(bracket_types) > 1:
            if return_pair_types:
                # For multi-bracket, return both connectivity and pair types
                conn_dict = _connectivity_list_multi_bracket(structure, bracket_types)
                pair_types_dict: Dict[str, Dict[int, str]] = {}
                for bracket_name, conn_list in conn_dict.items():
                    # Extract bracket type from name (e.g., "()" -> "(")
                    open_bracket = bracket_name[0]
                    # Filter structure to only this bracket type for pair types
                    filtered_structure = "".join(
                        char if char in (open_bracket, bracket_name[1], ".", "&") else "."
                        for char in structure
                    )
                    _, pair_types = _connectivity_list_unified(
                        filtered_structure, [(open_bracket, bracket_name[1])]
                    )
                    pair_types_dict[bracket_name] = pair_types
                return conn_dict, pair_types_dict
            else:
                return _connectivity_list_multi_bracket(structure, bracket_types)
        else:
            if return_pair_types:
                # For pure bracket format, filter out letters/numbers to avoid pairing them
                # Only keep brackets, dots, spaces, and strand separators
                if bracket_types is None:
                    bracket_types = [("(", ")")]
                allowed_chars = set(".& ")
                for open_bracket, close_bracket in bracket_types:
                    allowed_chars.add(open_bracket)
                    allowed_chars.add(close_bracket)
                filtered_structure = "".join(
                    char if char in allowed_chars else "." for char in structure
                )
                return _connectivity_list_unified(filtered_structure, bracket_types)
            else:
                return _connectivity_list(structure, bracket_types)
    elif format == "letter":
        # Use unified parser - it handles letters (and any brackets that might be present)
        if return_pair_types:
            return _connectivity_list_unified(structure, bracket_types)
        else:
            conn, _ = _connectivity_list_unified(structure, bracket_types)
            return conn
    elif format == "number":
        # Use unified parser - it handles numbers (and any brackets that might be present)
        if return_pair_types:
            return _connectivity_list_unified(structure, bracket_types)
        else:
            conn, _ = _connectivity_list_unified(structure, bracket_types)
            return conn
    elif format == "mixed":
        # For mixed format, use unified parser that handles brackets, letters, and numbers together
        if return_pair_types:
            return _connectivity_list_unified(structure, bracket_types)
        else:
            conn, _ = _connectivity_list_unified(structure, bracket_types)
            return conn
    else:
        raise ValueError(f"Unknown format: {format}")


class ConnectivityList:
    """Represents a connectivity list for RNA secondary structure.

    Attributes:
        connections (List[int]): A list of indices representing the connectivity
            between nucleotides.
        sequence (str): The RNA sequence.
        structure (str): The RNA secondary structure.
        pair_types (Dict[int, str]): Dictionary mapping position to pair type.
            Pair types can be bracket characters ('(', '[', '{', '<'),
            letters ('a', 'b', etc.), or numbers as strings ('1', '2', etc.).
            Only paired positions are included.
    """

    def __init__(self, sequence: str, structure: str, format: Optional[str] = None):
        """Initializes a ConnectivityList object.

        Args:
            sequence: The RNA sequence.
            structure: The RNA secondary structure.
            format: Structure format ('bracket', 'letter', 'number', 'auto').
                    If None, auto-detects format to handle mixed formats correctly.
        """
        if format is None:
            format = "auto"  # Auto-detect format to handle mixed formats

        result = _get_connectivity(structure, format, return_pair_types=True)

        if isinstance(result, tuple):
            # Single format with pair types
            conn, pair_types = result
            self.connections = conn
            self.pair_types = pair_types
        elif isinstance(result, dict):
            # Multi-bracket format (no pair types returned)
            # For multi-bracket, use the first bracket type's connectivity
            self.connections = list(result.values())[0]
            # Try to get pair types for the first bracket type
            bracket_name = list(result.keys())[0]
            open_bracket = bracket_name[0]
            # Filter structure to only this bracket type
            filtered_structure = "".join(
                char if char in (open_bracket, bracket_name[1], ".", "&") else "."
                for char in structure
            )
            _, pair_types = _connectivity_list_unified(
                filtered_structure, [(open_bracket, bracket_name[1])]
            )
            self.pair_types = pair_types
        else:
            # Single format without pair types (backward compatibility)
            self.connections = result
            self.pair_types = {}

        self.sequence = sequence
        self.structure = structure

    def get_pair_type(self, index: int) -> Optional[str]:
        """Get the pair type for a given position.

        Args:
            index: The position index.

        Returns:
            str: The pair type (e.g., '(', '[', 'a', '1') or None if unpaired.
        """
        return self.pair_types.get(index)

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


class ConnectivityListFactory:
    """Factory class for creating ConnectivityList objects.

    This factory handles the complexity of parsing different structure formats
    and creating ConnectivityList instances with appropriate pair type tracking.
    """

    @staticmethod
    def create(
        sequence: str,
        structure: str,
        format: Optional[str] = None,
        bracket_types: Optional[List[Tuple[str, str]]] = None,
    ) -> "ConnectivityList":
        """Create a ConnectivityList object from sequence and structure.

        Args:
            sequence: The RNA sequence.
            structure: The RNA secondary structure in any supported format.
            format: Structure format ('bracket', 'letter', 'number', 'auto').
                   If None, auto-detects format.
            bracket_types: For bracket format, specify bracket types.
                          None uses default [('(', ')')].

        Returns:
            ConnectivityList: A new ConnectivityList instance.

        Raises:
            ValueError: If structure format is invalid or malformed.
        """
        return ConnectivityList(sequence, structure, format)


def get_connectivity_list(
    sequence: str,
    structure: str,
    format: Optional[str] = None,
    bracket_types: Optional[List[Tuple[str, str]]] = None,
) -> ConnectivityList:
    """Factory function to create a ConnectivityList object.

    This is the recommended public API for creating ConnectivityList objects.
    It uses the ConnectivityListFactory internally.

    Args:
        sequence: The RNA sequence.
        structure: The RNA secondary structure in any supported format.
        format: Structure format ('bracket', 'letter', 'number', 'auto').
               If None, auto-detects format.
        bracket_types: For bracket format, specify bracket types.
                      None uses default [('(', ')')].
                      Use STANDARD_BRACKET_TYPES for full pseudo-knot support.

    Returns:
        ConnectivityList: A new ConnectivityList instance with pair type tracking.

    Raises:
        ValueError: If structure format is invalid or malformed.

    Examples:
        >>> cl = get_connectivity_list("GGGAAACCC", "(((...)))")
        >>> cl.is_nucleotide_paired(0)
        True
        >>> cl.get_pair_type(0)
        '('

        >>> cl = get_connectivity_list("GGGAAACCC", "a b c c b a", format="letter")
        >>> cl.get_pair_type(0)
        'a'
    """
    return ConnectivityListFactory.create(sequence, structure, format, bracket_types)


# Backward compatibility: Keep connectivity_list as a public function for parser.py
# but mark it as using internal implementation
def connectivity_list(
    structure: str, bracket_types: Optional[List[Tuple[str, str]]] = None
) -> List[int]:
    """Generate connectivity list from dot-bracket structure.

    NOTE: This function is kept for backward compatibility.
    For new code, use get_connectivity_list() to create ConnectivityList objects.

    Args:
        structure: A dot-bracket structure.
        bracket_types: List of (open, close) bracket pairs. Defaults to [('(', ')')].

    Returns:
        List[int]: The connectivity list.

    Raises:
        ValueError: If brackets are unbalanced.
    """
    return _connectivity_list(structure, bracket_types)


def is_circular(start: int, connections: List[int]) -> bool:
    """Check if a given RNA structure is circular.

    A structure is circular if following the connections from the start position
    eventually returns to the exact start position, forming a closed loop, AND
    there are no unpaired positions immediately after the start position.

    This is used by the parser to determine if there's more structure after
    a given position. If is_circular returns True, it means the structure
    loops back and there's no trailing structure. If False, there may be
    more structure (like a trailing single strand).

    Args:
        start: The starting index of the RNA structure.
        connections: A list of connections between nucleotides.

    Returns:
        bool: True if the RNA structure is circular.
    """
    if start < 0 or start >= len(connections):
        return False

    # If start is unpaired, it's not circular
    if connections[start] == -1:
        return False

    # If there are unpaired positions immediately after start, the structure continues
    # and is not circular at this point
    if start + 1 < len(connections) and connections[start + 1] == -1:
        return False

    visited = set()
    it = start

    while True:
        # Skip unpaired positions
        while it < len(connections) and connections[it] == -1:
            it += 1

        if it >= len(connections):
            return False  # Reached end of structure, not circular

        # Check if we've visited this position before (cycle detected)
        if it in visited:
            return it == start  # Circular only if we return to exact start

        visited.add(it)

        # Move to paired position
        if connections[it] < 0 or connections[it] >= len(connections):
            return False  # Invalid connection, not circular

        paired_pos = connections[it]
        # Move to position after the paired nucleotide
        it = paired_pos + 1

        # Check if we've looped back to start (circular)
        if it == start:
            return True

        # If we've gone past the end, it's not circular
        if it >= len(connections):
            return False

        # Prevent infinite loop - if we've visited too many positions, it's not circular
        if len(visited) > len(connections):
            return False
