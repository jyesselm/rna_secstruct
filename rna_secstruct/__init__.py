__author__ = "Joe Yesselman"
__email__ = "jyesselm@unl.edu"
__version__ = "0.1.0"

from .secstruct import SecStruct, MotifSearchParams
from .connectivity import (
    get_connectivity_list,
    ConnectivityList,
    STANDARD_BRACKET_TYPES,
)

# Auto-register pandas extensions if pandas is available
try:
    from . import pandas_extensions  # noqa: F401
except ImportError:
    pass  # pandas not available

# Auto-register parallel module
try:
    from . import parallel  # noqa: F401
except ImportError:
    pass  # parallel dependencies not available
