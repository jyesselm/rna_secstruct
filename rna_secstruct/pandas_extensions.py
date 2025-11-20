"""Pandas integration for RNA secondary structures.

This module provides pandas accessors (`.rna`) for working with SecStruct
objects in DataFrames and Series. No pandas modification is required - we
use pandas' extension API to register accessors.
"""

try:
    import pandas as pd
    import json
    from typing import Optional, List, Union
    from rna_secstruct.secstruct import SecStruct
    from rna_secstruct.json_encoder import SecStructJSONEncoder

    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False
    pd = None


if PANDAS_AVAILABLE:
    # Register custom encoder with pandas for JSON serialization
    try:
        # Try to patch pandas JSON encoder if available
        if hasattr(pd.io, "json") and hasattr(pd.io.json, "_json"):
            _original_default = getattr(pd.io.json._json.JSONEncoder, "default", None)

            def _patched_default(self, obj):
                """Patched default method for pandas JSON encoder."""
                if isinstance(obj, SecStruct):
                    return obj.to_dict()
                elif hasattr(obj, "to_dict") and callable(getattr(obj, "to_dict", None)):
                    # Check if it's a Motif (avoid circular import)
                    try:
                        from rna_secstruct.motif import Motif

                        if isinstance(obj, Motif):
                            return obj.to_dict()
                    except ImportError:
                        pass
                if _original_default:
                    return _original_default(self, obj)
                return super(pd.io.json._json.JSONEncoder, self).default(obj)

            # Only patch if we can access the encoder
            if _original_default:
                pd.io.json._json.JSONEncoder.default = _patched_default
    except (AttributeError, ImportError):
        pass  # pandas JSON encoder not available or different version

    @pd.api.extensions.register_dataframe_accessor("rna")
    class RNADataFrameAccessor:
        """Pandas DataFrame accessor for RNA structures.

        This is automatically available on all DataFrames as .rna
        after importing rna_secstruct. No pandas modification needed.
        """

        def __init__(self, pandas_obj):
            """Initialize accessor.

            Args:
                pandas_obj: The DataFrame that .rna was called on.
            """
            self._obj = pandas_obj

        def from_sequence_structure(
            self, seq_col: str, struct_col: str
        ) -> "pd.Series":
            """Create SecStruct objects from sequence and structure columns.

            Args:
                seq_col: Column name with sequences.
                struct_col: Column name with structures.

            Returns:
                pandas.Series: Series of SecStruct objects.
            """
            return self._obj.apply(
                lambda row: SecStruct(row[seq_col], row[struct_col]), axis=1
            )

        def add_secstruct(
            self, seq_col: str, struct_col: str, column: str = "secstruct"
        ) -> "pd.DataFrame":
            """Add SecStruct column to DataFrame.

            Args:
                seq_col: Column name with sequences.
                struct_col: Column name with structures.
                column: Name for the new SecStruct column.

            Returns:
                pd.DataFrame: DataFrame with added SecStruct column.
            """
            df = self._obj.copy()
            df[column] = self.from_sequence_structure(seq_col, struct_col)
            return df

        def add_statistics(self, secstruct_col: str = "secstruct") -> "pd.DataFrame":
            """Add multiple statistics columns to DataFrame.

            Args:
                secstruct_col: Column name with SecStruct objects.

            Returns:
                pd.DataFrame: DataFrame with added statistics columns.
            """
            df = self._obj.copy()
            df[f"{secstruct_col}_num_bp"] = df[secstruct_col].apply(
                lambda s: s.get_num_basepairs()
            )
            df[f"{secstruct_col}_num_unpaired"] = df[secstruct_col].apply(
                lambda s: s.get_num_unpaired()
            )
            df[f"{secstruct_col}_gc_content"] = df[secstruct_col].apply(
                lambda s: s.get_gc_content()
            )
            df[f"{secstruct_col}_length"] = df[secstruct_col].apply(len)
            return df

    @pd.api.extensions.register_series_accessor("rna")
    class RNASeriesAccessor:
        """Pandas Series accessor for RNA structures.

        This is automatically available on Series containing SecStruct objects.
        """

        def __init__(self, pandas_obj):
            """Initialize accessor.

            Args:
                pandas_obj: The Series that .rna was called on.
            """
            self._obj = pandas_obj

        def to_json(self, **kwargs) -> str:
            """Convert Series of SecStruct to JSON.

            Args:
                **kwargs: Arguments for json.dumps.

            Returns:
                str: JSON string representation.
            """
            structs = [s.to_dict() if isinstance(s, SecStruct) else s for s in self._obj]
            return json.dumps(structs, cls=SecStructJSONEncoder, **kwargs)

        def from_json(self, json_str: str) -> "pd.Series":
            """Create Series of SecStruct from JSON.

            Args:
                json_str: JSON string representation.

            Returns:
                pd.Series: Series of SecStruct objects.
            """
            data = json.loads(json_str)
            structs = [SecStruct.from_dict(d) for d in data]
            return pd.Series(structs, index=self._obj.index)

        def num_basepairs(self) -> "pd.Series":
            """Get number of base pairs for each structure.

            Returns:
                pd.Series: Series of base pair counts.
            """
            return self._obj.apply(
                lambda s: s.get_num_basepairs() if isinstance(s, SecStruct) else None
            )

        def num_motifs(self) -> "pd.Series":
            """Get number of motifs for each structure.

            Returns:
                pd.Series: Series of motif counts.
            """
            return self._obj.apply(
                lambda s: s.get_num_motifs() if isinstance(s, SecStruct) else None
            )

        def gc_content(self) -> "pd.Series":
            """Get GC content for each structure.

            Returns:
                pd.Series: Series of GC content values.
            """
            return self._obj.apply(
                lambda s: s.get_gc_content() if isinstance(s, SecStruct) else None
            )

        def helix_lengths(self) -> "pd.Series":
            """Get helix lengths for each structure.

            Returns:
                pd.Series: Series of helix length lists.
            """
            return self._obj.apply(
                lambda s: s.get_helix_lengths() if isinstance(s, SecStruct) else None
            )

        def has_pseudoknot(self) -> "pd.Series":
            """Check if structures have pseudo-knots.

            Returns:
                pd.Series: Series of boolean values.
            """
            return self._obj.apply(
                lambda s: (
                    hasattr(s, "connectivity")
                    and len([c for c in s.connectivity if c != -1]) > 0
                )
                if isinstance(s, SecStruct)
                else None
            )

