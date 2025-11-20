"""Tests for parallel processing."""

import pytest
from rna_secstruct.parallel import batch_parse, batch_connectivity, batch_apply
from rna_secstruct.secstruct import SecStruct


class TestParallelProcessing:
    """Test parallel processing functions."""

    def test_batch_parse_sequential(self):
        """Test batch parsing with sequential backend."""
        sequences = ["GGGAAACCC", "AAAGGGCCC", "CCCGAAAGGG"]
        structures = ["(((...)))", "(((...)))", "(((...)))"]  # All 9 chars - need to match
        # Fix: make all sequences 9 chars
        sequences = ["GGGAAACCC", "AAAGGGCCC", "CCCGAAAGG"]  # Last one is 9 chars
        result = batch_parse(sequences, structures, n_jobs=1, backend="sequential")
        assert len(result) == 3
        assert all(isinstance(s, SecStruct) for s in result)
        assert result[0].sequence == "GGGAAACCC"
        assert result[1].sequence == "AAAGGGCCC"
        assert result[2].sequence == "CCCGAAAGG"

    def test_batch_parse_length_mismatch(self):
        """Test batch_parse with mismatched lengths."""
        sequences = ["GGGAAACCC", "AAAGGGCCC"]
        structures = ["(((...)))"]
        with pytest.raises(ValueError, match="must have the same length"):
            batch_parse(sequences, structures)

    def test_batch_connectivity_sequential(self):
        """Test batch connectivity with sequential backend."""
        sequences = ["GGGAAACCC", "AAAGGGCCC"]
        structures = ["(((...)))", "(((...)))"]
        result = batch_connectivity(
            sequences, structures, n_jobs=1, backend="sequential"
        )
        assert len(result) == 2
        assert all(hasattr(r, "connections") or isinstance(r, list) for r in result)

    def test_batch_connectivity_length_mismatch(self):
        """Test batch_connectivity with mismatched lengths."""
        sequences = ["GGGAAACCC", "AAAGGGCCC"]
        structures = ["(((...)))"]
        with pytest.raises(ValueError, match="must have the same length"):
            batch_connectivity(sequences, structures)

    def test_batch_apply_sequential(self):
        """Test batch apply with sequential backend."""
        structs = [
            SecStruct("GGGAAACCC", "(((...)))"),
            SecStruct("AAAGGGCCC", "(((...)))"),
        ]
        result = batch_apply(
            structs, lambda s: s.get_num_basepairs(), n_jobs=1, backend="sequential"
        )
        assert len(result) == 2
        assert all(isinstance(n, int) for n in result)

    def test_batch_parse_invalid_backend(self):
        """Test batch_parse with invalid backend."""
        sequences = ["GGGAAACCC"]
        structures = ["(((...)))"]
        with pytest.raises(ValueError, match="Unknown backend"):
            batch_parse(sequences, structures, backend="invalid")

    def test_batch_parse_n_jobs_none(self):
        """Test batch_parse with n_jobs=None (auto-detect)."""
        sequences = ["GGGAAACCC", "AAAGGGCCC"]
        structures = ["(((...)))", "(((...)))"]
        result = batch_parse(sequences, structures, n_jobs=None, backend="sequential")
        assert len(result) == 2

