from enum import Enum

from .StandardWFOMC import standard_wfomc
from .FastWFOMC import fast_wfomc
from .IncrementalWFOMC import incremental_wfomc
from .RecursiveWFOMC import recursive_wfomc
from .TDWFOMC import td_wfomc

__all__ = [
    "standard_wfomc",
    "fast_wfomc",
    "incremental_wfomc",
    "recursive_wfomc",
    "tree"
]


class Algo(Enum):
    STANDARD = 'standard'
    FAST = 'fast'
    FASTv2 = 'fastv2'
    INCREMENTAL = 'incremental'
    RECURSIVE = 'recursive'
    TREE_DECOMP = "tree"

    def __str__(self):
        return self.value
