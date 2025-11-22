"""Custom JSON encoder for SecStruct and related objects."""

import json
from typing import Any

try:
    from rna_secstruct.secstruct import SecStruct
    from rna_secstruct.motif import Motif
except ImportError:
    SecStruct = None
    Motif = None


class SecStructJSONEncoder(json.JSONEncoder):
    """Custom JSON encoder for SecStruct and related objects.

    Automatically handles SecStruct and Motif objects in JSON serialization.
    """

    def default(self, obj: Any) -> Any:
        """Convert objects to JSON-serializable format.

        Args:
            obj: Object to encode.

        Returns:
            JSON-serializable representation of the object.
        """
        if SecStruct is not None and isinstance(obj, SecStruct):
            return obj.to_dict()
        elif Motif is not None and isinstance(obj, Motif):
            return obj.to_dict()
        # Let base class handle other types
        return super().default(obj)


def dumps(obj: Any, **kwargs) -> str:
    """JSON dumps with SecStruct support.

    Args:
        obj: Object to serialize.
        **kwargs: Additional arguments for json.dumps.

    Returns:
        str: JSON string representation.
    """
    return json.dumps(obj, cls=SecStructJSONEncoder, **kwargs)


def loads(json_str: str, **kwargs) -> Any:
    """JSON loads (standard, but provided for symmetry).

    Args:
        json_str: JSON string to parse.
        **kwargs: Additional arguments for json.loads.

    Returns:
        Parsed JSON object.
    """
    return json.loads(json_str, **kwargs)
