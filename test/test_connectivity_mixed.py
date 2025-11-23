"""
Tests for mixed format connectivity (brackets + letters/numbers).
"""

from rna_secstruct.connectivity import (
    _detect_structure_format,
    get_connectivity_list,
)


class TestMixedFormats:
    """Test mixed format structures with brackets and letters/numbers."""

    def test_mixed_brackets_and_letters(self):
        """Test structure with brackets and letters mixed together."""
        # Structure and sequence must have same length
        structure = "(((aaa(((...)))aaa)))"
        seq = "GGGAAACCCUUUAAAGGGCCC"  # 21 nucleotides to match structure
        cl = get_connectivity_list(seq, structure)

        # Should detect as mixed format
        assert _detect_structure_format(structure) == "mixed"

        # Check bracket pairs
        assert cl.get_pair_type(0) == "("  # Opening bracket
        assert cl.get_pair_type(20) == ")"  # Closing bracket
        assert cl.get_pair_type(6) == "("  # Inner opening bracket
        assert cl.get_pair_type(14) == ")"  # Inner closing bracket

        # Check letter pairs
        assert cl.get_pair_type(3) == "a"  # First 'a'
        assert cl.get_pair_type(17) == "a"  # Paired 'a'
        assert cl.get_pair_type(4) == "a"  # Second 'a'
        assert cl.get_pair_type(16) == "a"  # Paired 'a'
        assert cl.get_pair_type(5) == "a"  # Third 'a'
        assert cl.get_pair_type(15) == "a"  # Paired 'a'

        # Verify connections
        assert cl.connections[3] == 17  # 'a' at 3 pairs with 'a' at 17
        assert cl.connections[4] == 16  # 'a' at 4 pairs with 'a' at 16
        assert cl.connections[5] == 15  # 'a' at 5 pairs with 'a' at 15

    def test_letters_with_brackets(self):
        """Test structure with letters and brackets."""
        structure = "aaa(((...)))aaa"
        seq = "GGGAAACCCUUUAAA"  # 15 nucleotides to match structure
        cl = get_connectivity_list(seq, structure)

        # Should detect as mixed format
        assert _detect_structure_format(structure) == "mixed"

        # Check bracket pairs
        assert cl.get_pair_type(3) == "("  # Opening bracket
        assert cl.get_pair_type(9) == ")"  # Closing bracket

        # Check letter pairs
        assert cl.get_pair_type(0) == "a"
        assert cl.get_pair_type(12) == "a"
        assert cl.get_pair_type(1) == "a"
        assert cl.get_pair_type(13) == "a"
        assert cl.get_pair_type(2) == "a"
        assert cl.get_pair_type(14) == "a"

    def test_mixed_brackets_and_numbers(self):
        """Test structure with brackets and numbers."""
        structure = "(((111(((...)))111)))"
        seq = "GGGAAACCCUUUAAAGGGCCC"  # 21 nucleotides to match structure
        cl = get_connectivity_list(seq, structure)

        # Should detect as mixed format
        assert _detect_structure_format(structure) == "mixed"

        # Check bracket pairs exist
        assert cl.get_pair_type(0) == "("
        assert cl.get_pair_type(20) == ")"

        # Check number pairs
        # Numbers are at positions where digits start
        # '1' at position 3, 4, 5 should pair with '1' at 15, 16, 17
        # But numbers are multi-digit, so we need to check the actual positions
        num_pair_types = {k: v for k, v in cl.pair_types.items() if v.isdigit()}
        assert len(num_pair_types) > 0  # Should have number pair types

    def test_pure_bracket_with_letters_ignored(self):
        """Test that pure bracket format ignores letters when not in mixed mode."""
        structure = "((x...))"
        seq = "GGGAAACC"  # 8 nucleotides to match structure
        cl = get_connectivity_list(seq, structure, format="bracket")

        # 'x' should be ignored (not treated as a letter to pair)
        assert cl.get_pair_type(2) is None  # 'x' is unpaired
        assert cl.connections[2] == -1  # 'x' is unpaired
        # Brackets should still work
        assert cl.get_pair_type(0) == "("
        assert cl.get_pair_type(7) == ")"
