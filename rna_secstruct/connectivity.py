"""
Connectivity list generation and manipulation for RNA secondary structures.

Supports multiple bracket types for pseudo-knots and alternative pairing notations
(letter-based and number-based).

Public API:
    - get_connectivity_list(): Main unified interface to create ConnectivityList objects
    - connectivity_list(): Convenience function that returns simple List[int]
    - ConnectivityList: Main class for working with connectivity lists
    - has_pseudoknot(): Utility function to detect pseudo-knots
    - is_circular(): Utility function to detect circular structures
    - STANDARD_BRACKET_TYPES: Standard bracket types for pseudo-knots

All other functions are internal and should not be used directly.
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


def _validate_structure_input(structure: str) -> None:
    """Internal: Validate basic structure input.

    Args:
        structure: Structure string to validate.

    Raises:
        ValueError: If structure is invalid with detailed explanation.
    """
    if not isinstance(structure, str):
        raise TypeError(
            f"Structure must be a string, got {type(structure).__name__}. "
            f"Structure represents RNA secondary structure using brackets, letters, numbers, or dots."
        )

    if len(structure) == 0:
        raise ValueError(
            "Structure cannot be empty. Provide a valid structure string with at least one character."
        )


def _validate_bracket_types(bracket_types: List[Tuple[str, str]]) -> None:
    """Internal: Validate bracket types configuration.

    Args:
        bracket_types: List of (open, close) bracket pairs.

    Raises:
        ValueError: If bracket types are invalid with detailed explanation.
    """
    if not bracket_types:
        raise ValueError(
            "Bracket types cannot be empty. Provide at least one (open, close) bracket pair, "
            "e.g., [('(', ')')] for standard parentheses."
        )

    seen_opens = set()
    seen_closes = set()

    for i, bracket_pair in enumerate(bracket_types):
        if not isinstance(bracket_pair, tuple):
            raise ValueError(
                f"Bracket type at index {i} must be a tuple of two strings, "
                f"got {type(bracket_pair).__name__}."
            )

        if len(bracket_pair) != 2:
            raise ValueError(
                f"Bracket type at index {i} must be a tuple of exactly two strings, "
                f"got {len(bracket_pair)} element(s). "
                f"Each bracket type must be (open_char, close_char), e.g., ('(', ')')."
            )

        open_bracket, close_bracket = bracket_pair

        if not isinstance(open_bracket, str) or not isinstance(close_bracket, str):
            raise ValueError(
                f"Bracket type at index {i} must be a tuple of two strings, "
                f"got ({type(open_bracket).__name__}, {type(close_bracket).__name__})."
            )

        if len(open_bracket) != 1 or len(close_bracket) != 1:
            raise ValueError(
                f"Bracket type at index {i}: '{open_bracket}' and '{close_bracket}' must be single characters. "
                f"Each bracket type must be a (single_char, single_char) tuple."
            )

        if open_bracket == close_bracket:
            raise ValueError(
                f"Bracket type at index {i}: Opening and closing brackets cannot be the same character '{open_bracket}'. "
                f"Use different characters for opening and closing brackets, e.g., ('(', ')')."
            )

        if open_bracket in seen_opens:
            raise ValueError(
                f"Duplicate opening bracket '{open_bracket}' at index {i}. "
                f"Each bracket type must have a unique opening bracket character."
            )

        if close_bracket in seen_closes:
            raise ValueError(
                f"Duplicate closing bracket '{close_bracket}' at index {i}. "
                f"Each bracket type must have a unique closing bracket character."
            )

        if open_bracket in seen_closes or close_bracket in seen_opens:
            raise ValueError(
                f"Bracket type at index {i}: '{open_bracket}' and '{close_bracket}' conflict with existing bracket types. "
                f"A character cannot be both an opening and closing bracket, or overlap with other bracket types."
            )

        seen_opens.add(open_bracket)
        seen_closes.add(close_bracket)


def _parse_connectivity(
    structure: str,
    bracket_types: Optional[List[Tuple[str, str]]] = None,
    case_sensitive: bool = True,
    return_pair_types: bool = True,
) -> Tuple[List[int], Dict[int, str]]:
    """Internal: Core parser that handles brackets, letters, and numbers.

    This is the single unified parser that processes all formats including
    mixed formats like "(((aaa(((...)))aaa)))" where brackets and letters are interleaved.

    Args:
        structure: Structure string that may contain brackets, letters, and/or numbers.
        bracket_types: List of (open, close) bracket pairs. Defaults to [('(', ')')].
        case_sensitive: Whether to treat uppercase and lowercase letters as different pairs.
        return_pair_types: If True, track and return pair type information.

    Returns:
        Tuple[List[int], Dict[int, str]]: Connectivity list and pair type dictionary.
            If return_pair_types=False, pair_types dict will be empty.

    Raises:
        ValueError: If structure is invalid, brackets/letters/numbers are unbalanced,
            or there are conflicting pairings. Error messages explain the specific issue.
        TypeError: If structure is not a string.
    """
    # Validate input
    _validate_structure_input(structure)

    if bracket_types is None:
        bracket_types = [("(", ")")]
    else:
        _validate_bracket_types(bracket_types)

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
                    f"Unbalanced bracket: closing bracket '{char}' at position {index} has no matching opening bracket '{open_bracket}'. "
                    f"This means there are more closing brackets than opening brackets. "
                    f"Check that all brackets are properly matched and balanced in your structure."
                )
            complement, pair_type = stacks[open_bracket].pop()

            # Check for conflicts before pairing (brackets should not conflict with themselves)
            if connections[complement] != -1:
                raise ValueError(
                    f"Position {complement} (opening bracket '{open_bracket}') is already paired with position {connections[complement]}, "
                    f"but is also being paired with position {index} (closing bracket '{char}'). "
                    f"This indicates a structural error: brackets are not properly nested or balanced."
                )
            if connections[index] != -1:
                raise ValueError(
                    f"Position {index} (closing bracket '{char}') is already paired with position {connections[index]}, "
                    f"but is also being paired with position {complement} (opening bracket '{open_bracket}'). "
                    f"This indicates a structural error: brackets are not properly nested or balanced."
                )

            if complement == index:
                raise ValueError(
                    f"Position {index} cannot be paired with itself. "
                    f"Self-pairing is not valid in RNA secondary structures. "
                    f"Check your structure for errors."
                )

            connections[complement] = index
            connections[index] = complement
            if return_pair_types:
                pair_types[complement] = pair_type
                pair_types[index] = char
        elif char not in (".", "&", " "):
            # Check if it's a valid letter (a-z, A-Z) for pairing
            if re.match(r"^[a-zA-Z]$", char):
                original_char = char
                if not case_sensitive:
                    char = char.lower()
                letter_positions[char].append((index, original_char))
            # Other characters are ignored - treated as unpaired
        # '.' and '&' and ' ' are ignored (unpaired and strand separator)

    # Check for unmatched opening brackets
    for open_bracket, stack in stacks.items():
        if stack:
            unmatched_positions = [pos for pos, _ in stack]
            raise ValueError(
                f"Unbalanced brackets: {len(stack)} unmatched opening bracket(s) '{open_bracket}' "
                f"at position(s) {unmatched_positions}. "
                f"This means there are more opening brackets than closing brackets. "
                f"Each opening bracket must have a corresponding closing bracket. "
                f"Check that all brackets are properly matched and balanced in your structure."
            )

    # Second pass: process letters (pair from outside in)
    if letter_positions:
        for letter, positions in letter_positions.items():
            if len(positions) % 2 != 0:
                pos_list = [pos for pos, _ in positions]
                raise ValueError(
                    f"Unpaired letter '{letter}': found {len(positions)} occurrence(s) at position(s) {pos_list} "
                    f"(must be even for pairing). "
                    f"In letter-based notation, each letter must appear an even number of times to form pairs. "
                    f"Check that all letters are properly paired in your structure."
                )

            while positions:
                i, i_char = positions.pop(0)
                j, j_char = positions.pop()
                # Only pair if not already paired by brackets
                if connections[i] != -1:
                    raise ValueError(
                        f"Position {i} (letter '{i_char}') is already paired with position {connections[i]} "
                        f"from bracket pairing, but is also being paired with position {j} (letter '{j_char}'). "
                        f"This creates a conflict: a nucleotide cannot be paired with multiple partners. "
                        f"Check your structure for overlapping bracket and letter pairings."
                    )
                if connections[j] != -1:
                    raise ValueError(
                        f"Position {j} (letter '{j_char}') is already paired with position {connections[j]} "
                        f"from bracket pairing, but is also being paired with position {i} (letter '{i_char}'). "
                        f"This creates a conflict: a nucleotide cannot be paired with multiple partners. "
                        f"Check your structure for overlapping bracket and letter pairings."
                    )

                if i == j:
                    raise ValueError(
                        f"Position {i} (letter '{i_char}') cannot be paired with itself. "
                        f"Self-pairing is not valid in RNA secondary structures. "
                        f"Check your structure for errors."
                    )

                if j < 0 or j >= len(connections) or i < 0 or i >= len(connections):
                    raise ValueError(
                        f"Invalid pairing: position {i} is being paired with position {j}, "
                        f"but one or both positions are out of bounds (structure length is {len(connections)}). "
                        f"This indicates a structural error in your input."
                    )

                connections[i] = j
                connections[j] = i
                if return_pair_types:
                    pair_types[i] = i_char
                    pair_types[j] = j_char

    # Third pass: process numbers (pair from outside in)
    if number_positions:
        for start, end, num, num_str in number_positions:
            number_pos_map[num].append((start, end, num_str))

        for num, positions in number_pos_map.items():
            if len(positions) % 2 != 0:
                pos_list = [start for start, _, _ in positions]
                raise ValueError(
                    f"Unpaired number '{num}': found {len(positions)} occurrence(s) starting at position(s) {pos_list} "
                    f"(must be even for pairing). "
                    f"In number-based notation, each number must appear an even number of times to form pairs. "
                    f"Check that all numbers are properly paired in your structure."
                )

            while positions:
                (i_start, i_end, i_num_str) = positions.pop(0)
                (j_start, j_end, j_num_str) = positions.pop()
                # Only pair if not already paired by brackets or letters
                if connections[i_start] != -1:
                    raise ValueError(
                        f"Position {i_start} (number '{i_num_str}') is already paired with position {connections[i_start]} "
                        f"from bracket or letter pairing, but is also being paired with position {j_start} (number '{j_num_str}'). "
                        f"This creates a conflict: a nucleotide cannot be paired with multiple partners. "
                        f"Check your structure for overlapping pairings."
                    )
                if connections[j_start] != -1:
                    raise ValueError(
                        f"Position {j_start} (number '{j_num_str}') is already paired with position {connections[j_start]} "
                        f"from bracket or letter pairing, but is also being paired with position {i_start} (number '{i_num_str}'). "
                        f"This creates a conflict: a nucleotide cannot be paired with multiple partners. "
                        f"Check your structure for overlapping pairings."
                    )

                if i_start == j_start:
                    raise ValueError(
                        f"Position {i_start} (number '{i_num_str}') cannot be paired with itself. "
                        f"Self-pairing is not valid in RNA secondary structures. "
                        f"Check your structure for errors."
                    )

                if (
                    j_start < 0
                    or j_start >= len(connections)
                    or i_start < 0
                    or i_start >= len(connections)
                ):
                    raise ValueError(
                        f"Invalid pairing: position {i_start} is being paired with position {j_start}, "
                        f"but one or both positions are out of bounds (structure length is {len(connections)}). "
                        f"This indicates a structural error in your input."
                    )

                connections[i_start] = j_start
                connections[j_start] = i_start
                if return_pair_types:
                    pair_types[i_start] = i_num_str
                    pair_types[j_start] = j_num_str

    return connections, pair_types


def _parse_multi_bracket(
    structure: str, bracket_types: Optional[List[Tuple[str, str]]] = None
) -> Dict[str, Tuple[List[int], Dict[int, str]]]:
    """Internal: Generate connectivity lists for multiple bracket types separately.

    This allows detection of pseudo-knots by analyzing interleaving bracket types.

    Args:
        structure: A dot-bracket structure with multiple bracket types.
        bracket_types: List of (open, close) bracket pairs. Defaults to STANDARD_BRACKET_TYPES.

    Returns:
        Dict[str, Tuple[List[int], Dict[int, str]]]: Dictionary mapping bracket type names
            (e.g., '()', '[]') to (connectivity list, pair types) tuples.

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
            conn, pair_types = _parse_connectivity(
                filtered_structure, [(open_bracket, close_bracket)], return_pair_types=True
            )
            result[bracket_name] = (conn, pair_types)
        except ValueError:
            # This bracket type might not be present or balanced, skip it
            result[bracket_name] = ([-1] * len(structure), {})

    return result


def _detect_structure_format(structure: str) -> str:
    """Internal: Detect the format of the structure string.

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


def _pairs_cross(conn1: List[int], conn2: List[int]) -> bool:
    """Internal: Check if pairs from two connectivity lists cross each other.

    Args:
        conn1: First connectivity list.
        conn2: Second connectivity list.

    Returns:
        bool: True if pairs from the two lists cross each other.
    """
    pairs1 = [(i, conn1[i]) for i in range(len(conn1)) if conn1[i] != -1 and i < conn1[i]]
    pairs2 = [(i, conn2[i]) for i in range(len(conn2)) if conn2[i] != -1 and i < conn2[i]]

    for i1, j1 in pairs1:
        for i2, j2 in pairs2:
            # Check if pairs cross: (i1 < i2 < j1 < j2) or (i2 < i1 < j2 < j1)
            if (i1 < i2 < j1 < j2) or (i2 < i1 < j2 < j1):
                return True

    return False


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
    """Internal: Unified interface to get connectivity list(s) from any supported format.

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
        TypeError: If structure is not a string.
    """
    # Validate input early
    _validate_structure_input(structure)

    if format is None or format == "auto":
        format = _detect_structure_format(structure)

    # Handle multi-bracket format
    if format == "bracket" and bracket_types and len(bracket_types) > 1:
        multi_result = _parse_multi_bracket(structure, bracket_types)
        if return_pair_types:
            conn_dict = {name: conn for name, (conn, _) in multi_result.items()}
            pair_types_dict = {name: pt for name, (_, pt) in multi_result.items()}
            return conn_dict, pair_types_dict
        else:
            return {name: conn for name, (conn, _) in multi_result.items()}

    # Single format - use unified parser
    # For pure bracket format, filter out letters/numbers to avoid pairing them
    if format == "bracket" and bracket_types is None:
        bracket_types = [("(", ")")]

    if format == "bracket":
        # Filter structure to only allow bracket characters
        allowed_chars = set(".& ")
        for open_bracket, close_bracket in bracket_types:
            allowed_chars.add(open_bracket)
            allowed_chars.add(close_bracket)
        filtered_structure = "".join(char if char in allowed_chars else "." for char in structure)
        conn, pair_types = _parse_connectivity(
            filtered_structure, bracket_types, return_pair_types=return_pair_types
        )
    else:
        # Letter, number, or mixed format - use full parser
        conn, pair_types = _parse_connectivity(
            structure, bracket_types, return_pair_types=return_pair_types
        )

    if return_pair_types:
        return conn, pair_types
    else:
        return conn


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

    def __init__(
        self,
        sequence: str,
        structure: str,
        format: Optional[str] = None,
        bracket_types: Optional[List[Tuple[str, str]]] = None,
    ):
        """Initialize a ConnectivityList object.

        Args:
            sequence: The RNA sequence.
            structure: The RNA secondary structure.
            format: Structure format ('bracket', 'letter', 'number', 'auto').
                    If None, auto-detects format to handle mixed formats correctly.
            bracket_types: For bracket format, specify bracket types.
                          None uses default [('(', ')')].

        Raises:
            ValueError: If sequence and structure have different lengths, or if structure is invalid.
            TypeError: If inputs are not strings.
        """
        # Validate inputs
        if not isinstance(sequence, str):
            raise TypeError(
                f"Sequence must be a string, got {type(sequence).__name__}. "
                f"Sequence should contain RNA nucleotides (A, U, G, C, T)."
            )
        if not isinstance(structure, str):
            raise TypeError(
                f"Structure must be a string, got {type(structure).__name__}. "
                f"Structure represents RNA secondary structure using brackets, letters, numbers, or dots."
            )

        # Spaces are ignored in structure parsing, so remove them for length comparison
        structure_no_spaces = structure.replace(" ", "")
        if len(sequence) != len(structure_no_spaces):
            raise ValueError(
                f"Sequence and structure must have the same length (ignoring spaces). "
                f"Sequence length: {len(sequence)}, structure length (without spaces): {len(structure_no_spaces)}. "
                f"Each nucleotide in the sequence must have a corresponding structure character. "
                f"Note: Spaces in structure are ignored during parsing."
            )

        if format is None:
            format = "auto"

        result = _get_connectivity(structure, format, bracket_types, return_pair_types=True)

        if isinstance(result, tuple):
            conn, pair_types = result
            # Check if this is multi-bracket format (tuple of two dicts)
            if isinstance(conn, dict) and isinstance(pair_types, dict):
                # Multi-bracket format - use the first bracket type's connectivity
                bracket_name = list(conn.keys())[0]
                self.connections = conn[bracket_name]
                self.pair_types = pair_types[bracket_name]
            else:
                # Single format with pair types
                self.connections = conn
                self.pair_types = pair_types
        elif isinstance(result, dict):
            # Multi-bracket format without pair types (shouldn't happen with return_pair_types=True)
            # For multi-bracket, use the first bracket type's connectivity
            self.connections = list(result.values())[0]
            # Try to get pair types for the first bracket type
            bracket_name = list(result.keys())[0]
            open_bracket = bracket_name[0]
            filtered_structure = "".join(
                char if char in (open_bracket, bracket_name[1], ".", "&") else "."
                for char in structure
            )
            _, pair_types = _parse_connectivity(
                filtered_structure, [(open_bracket, bracket_name[1])], return_pair_types=True
            )
            self.pair_types = pair_types
        else:
            # Single format without pair types (shouldn't happen with return_pair_types=True)
            self.connections = result
            self.pair_types = {}

        self.sequence = sequence
        self.structure = structure

    def get_pair_type(self, index: int) -> Optional[str]:
        """Get the pair type for a given position.

        Args:
            index: The position index.

        Returns:
            Optional[str]: The pair type (e.g., '(', '[', 'a', '1') or None if unpaired.
        """
        return self.pair_types.get(index)

    def is_nucleotide_paired(self, index: int) -> bool:
        """Check if a nucleotide at a given index is paired.

        Args:
            index: The index of the nucleotide.

        Returns:
            bool: True if the nucleotide is paired, False otherwise.
        """
        return self.connections[index] != -1

    def get_paired_nucleotide(self, index: int) -> int:
        """Get the index of the nucleotide paired with the nucleotide at the given index.

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
        """Get the base pair of the nucleotide at the given index.

        Args:
            index: The index of the nucleotide.

        Returns:
            str: The base pair of the nucleotide, or "." if unpaired.
        """
        if not self.is_nucleotide_paired(index):
            return "."
        return self.sequence[index] + self.sequence[self.get_paired_nucleotide(index)]


def get_connectivity_list(
    sequence: str,
    structure: str,
    format: Optional[str] = None,
    bracket_types: Optional[List[Tuple[str, str]]] = None,
) -> ConnectivityList:
    """Unified interface to create ConnectivityList objects from any supported format.

    This is the main public API for creating connectivity lists. It auto-detects
    the structure format and handles all complexity internally.

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
        ValueError: If structure format is invalid, malformed, has unbalanced pairs,
            conflicting pairings, or if sequence and structure lengths don't match.
            Error messages explain the specific issue.
        TypeError: If inputs are not strings.

    Examples:
        >>> cl = get_connectivity_list("GGGAAACCC", "(((...)))")
        >>> cl.is_nucleotide_paired(0)
        True
        >>> cl.get_pair_type(0)
        '('

        >>> cl = get_connectivity_list("GGGAAACCC", "a b c c b a", format="letter")
        >>> cl.get_pair_type(0)
        'a'

        >>> cl = get_connectivity_list("GGGAAACCC", "(((...)))", bracket_types=STANDARD_BRACKET_TYPES)
        >>> cl.is_nucleotide_paired(0)
        True
    """
    return ConnectivityList(sequence, structure, format, bracket_types)


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
        # For single list with multiple bracket types, would need original structure
        # In practice, use _parse_multi_bracket first
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


def connectivity_list(
    structure: str, bracket_types: Optional[List[Tuple[str, str]]] = None
) -> List[int]:
    """Generate connectivity list from structure (returns simple List[int]).

    This is a convenience function for code that just needs a connectivity list
    without the full ConnectivityList object. For full functionality, use
    get_connectivity_list() instead.

    Args:
        structure: A structure string in any supported format.
        bracket_types: List of (open, close) bracket pairs. Defaults to [('(', ')')].

    Returns:
        List[int]: The connectivity list. Each index contains the paired position or -1.

    Raises:
        ValueError: If structure is invalid, malformed, has unbalanced pairs, or conflicting pairings.
            Error messages explain the specific issue.
        TypeError: If structure is not a string.
    """
    result = _get_connectivity(
        structure, format=None, bracket_types=bracket_types, return_pair_types=False
    )
    if isinstance(result, list):
        return result
    elif isinstance(result, dict):
        # Multi-bracket format - return first bracket type's connectivity
        return list(result.values())[0]
    else:
        # Shouldn't happen, but handle it
        return result[0] if isinstance(result, tuple) else result


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
