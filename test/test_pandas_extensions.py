"""Tests for pandas integration."""

import pytest

try:
    import pandas as pd

    # Import to register accessors
    import rna_secstruct  # noqa: F401
    from rna_secstruct.secstruct import SecStruct

    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False
    pd = None


@pytest.mark.skipif(not PANDAS_AVAILABLE, reason="pandas not available")
class TestPandasExtensions:
    """Test pandas integration."""

    def test_dataframe_accessor_registered(self):
        """Test that DataFrame accessor is registered."""
        df = pd.DataFrame({"a": [1, 2]})
        assert hasattr(df, "rna")

    def test_series_accessor_registered(self):
        """Test that Series accessor is registered."""
        s = pd.Series([1, 2])
        assert hasattr(s, "rna")

    def test_from_sequence_structure(self):
        """Test creating SecStruct from DataFrame columns."""
        df = pd.DataFrame(
            {
                "sequence": ["GGGAAACCC", "AAAGGGCCC"],
                "structure": ["(((...)))", "(((...)))"],
            }
        )
        secstructs = df.rna.from_sequence_structure("sequence", "structure")
        assert len(secstructs) == 2
        assert all(isinstance(s, SecStruct) for s in secstructs)
        assert secstructs.iloc[0].sequence == "GGGAAACCC"

    def test_add_secstruct(self):
        """Test adding SecStruct column to DataFrame."""
        df = pd.DataFrame(
            {
                "sequence": ["GGGAAACCC", "AAAGGGCCC"],
                "structure": ["(((...)))", "(((...)))"],
            }
        )
        df_new = df.rna.add_secstruct("sequence", "structure")
        assert "secstruct" in df_new.columns
        assert all(isinstance(s, SecStruct) for s in df_new["secstruct"])

    def test_add_statistics(self):
        """Test adding statistics columns."""
        df = pd.DataFrame(
            {
                "sequence": ["GGGAAACCC", "AAAGGGCCC"],
                "structure": ["(((...)))", "(((...)))"],
            }
        )
        df["secstruct"] = df.rna.from_sequence_structure("sequence", "structure")
        df_new = df.rna.add_statistics("secstruct")
        assert "secstruct_num_bp" in df_new.columns
        assert "secstruct_num_unpaired" in df_new.columns
        assert "secstruct_gc_content" in df_new.columns
        assert "secstruct_length" in df_new.columns

    def test_series_num_basepairs(self):
        """Test Series accessor num_basepairs method."""
        structs = pd.Series(
            [
                SecStruct("GGGAAACCC", "(((...)))"),
                SecStruct("AAAGGGCCC", "(((...)))"),
            ]
        )
        num_bp = structs.rna.num_basepairs()
        assert len(num_bp) == 2
        assert all(isinstance(n, (int, type(None))) for n in num_bp)

    def test_series_gc_content(self):
        """Test Series accessor gc_content method."""
        structs = pd.Series(
            [
                SecStruct("GGGAAACCC", "(((...)))"),
                SecStruct("AAAGGGCCC", "(((...)))"),
            ]
        )
        gc = structs.rna.gc_content()
        assert len(gc) == 2
        assert all(isinstance(g, (float, type(None))) for g in gc)

    def test_series_to_json(self):
        """Test Series accessor to_json method."""
        structs = pd.Series(
            [
                SecStruct("GGGAAACCC", "(((...)))"),
                SecStruct("AAAGGGCCC", "(((...)))"),
            ]
        )
        json_str = structs.rna.to_json()
        assert isinstance(json_str, str)
        assert "GGGAAACCC" in json_str

    def test_series_from_json(self):
        """Test Series accessor from_json method."""
        structs = pd.Series(
            [
                SecStruct("GGGAAACCC", "(((...)))"),
                SecStruct("AAAGGGCCC", "(((...)))"),
            ]
        )
        json_str = structs.rna.to_json()
        structs2 = structs.rna.from_json(json_str)
        assert len(structs2) == 2
        assert all(isinstance(s, SecStruct) for s in structs2)
