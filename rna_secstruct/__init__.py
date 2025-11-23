__author__ = "Joe Yesselman"
__email__ = "jyesselm@unl.edu"
__version__ = "0.1.1"

import contextlib

from .connectivity import (
    STANDARD_BRACKET_TYPES,
    ConnectivityList,
    get_connectivity_list,
)
from .secstruct import MotifSearchParams, SecStruct

# Auto-register pandas extensions if pandas is available
with contextlib.suppress(ImportError):
    from . import pandas_extensions  # noqa: F401

# Auto-register parallel module
with contextlib.suppress(ImportError):
    from . import parallel  # noqa: F401
