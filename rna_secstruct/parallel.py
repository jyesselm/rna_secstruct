"""Parallel processing support for RNA secondary structures.

This module provides parallel batch operations for parsing and analyzing
multiple RNA structures efficiently.
"""

from typing import List, Optional, Callable, Any, Tuple
import multiprocessing as mp
from functools import partial

try:
    from rna_secstruct.secstruct import SecStruct
    from rna_secstruct.connectivity import get_connectivity_list
except ImportError:
    SecStruct = None
    get_connectivity_list = None


def _parse_single(sequence: str, structure: str) -> SecStruct:
    """Parse a single structure (helper for parallel processing).

    Args:
        sequence: RNA sequence.
        structure: RNA secondary structure.

    Returns:
        SecStruct: Parsed structure.
    """
    return SecStruct(sequence, structure)


def _get_connectivity_single(sequence: str, structure: str, format: Optional[str] = None) -> Any:
    """Get connectivity list for a single structure (helper for parallel processing).

    Args:
        sequence: RNA sequence.
        structure: RNA secondary structure.
        format: Structure format (optional).

    Returns:
        ConnectivityList or List[int]: Connectivity list.
    """
    return get_connectivity_list(sequence, structure, format=format)


def batch_parse(
    sequences: List[str],
    structures: List[str],
    n_jobs: Optional[int] = None,
    backend: str = "multiprocessing",
    **kwargs,
) -> List[SecStruct]:
    """Parse multiple structures in parallel.

    Args:
        sequences: List of RNA sequences.
        structures: List of RNA secondary structures.
        n_jobs: Number of parallel jobs. If None, uses CPU count.
        backend: Backend to use ('multiprocessing', 'threading', or 'sequential').
        **kwargs: Additional arguments passed to SecStruct constructor.

    Returns:
        List[SecStruct]: List of parsed structures.

    Raises:
        ValueError: If sequences and structures have different lengths.
    """
    if len(sequences) != len(structures):
        raise ValueError(
            f"Sequences and structures must have the same length. "
            f"Got {len(sequences)} sequences and {len(structures)} structures."
        )

    if n_jobs is None:
        n_jobs = mp.cpu_count()

    if backend == "sequential" or n_jobs == 1:
        # Sequential processing
        return [SecStruct(seq, struct, **kwargs) for seq, struct in zip(sequences, structures)]
    elif backend == "threading":
        # Threading backend (good for I/O-bound tasks, but GIL limits CPU-bound)
        from concurrent.futures import ThreadPoolExecutor

        with ThreadPoolExecutor(max_workers=n_jobs) as executor:
            return list(
                executor.map(
                    lambda args: SecStruct(args[0], args[1], **kwargs),
                    zip(sequences, structures),
                )
            )
    elif backend == "multiprocessing":
        # Multiprocessing backend (good for CPU-bound tasks)
        with mp.Pool(processes=n_jobs) as pool:
            return pool.starmap(
                lambda seq, struct: SecStruct(seq, struct, **kwargs),
                zip(sequences, structures),
            )
    else:
        raise ValueError(
            f"Unknown backend: {backend}. Must be 'multiprocessing', 'threading', or 'sequential'."
        )


def batch_connectivity(
    sequences: List[str],
    structures: List[str],
    format: Optional[str] = None,
    n_jobs: Optional[int] = None,
    backend: str = "multiprocessing",
    **kwargs,
) -> List[Any]:
    """Generate connectivity lists for multiple structures in parallel.

    Args:
        sequences: List of RNA sequences.
        structures: List of RNA secondary structures.
        format: Structure format (optional).
        n_jobs: Number of parallel jobs. If None, uses CPU count.
        backend: Backend to use ('multiprocessing', 'threading', or 'sequential').
        **kwargs: Additional arguments passed to get_connectivity_list.

    Returns:
        List: List of connectivity lists or ConnectivityList objects.

    Raises:
        ValueError: If sequences and structures have different lengths.
    """
    if len(sequences) != len(structures):
        raise ValueError(
            f"Sequences and structures must have the same length. "
            f"Got {len(sequences)} sequences and {len(structures)} structures."
        )

    if n_jobs is None:
        n_jobs = mp.cpu_count()

    if backend == "sequential" or n_jobs == 1:
        # Sequential processing
        return [
            get_connectivity_list(seq, struct, format=format, **kwargs)
            for seq, struct in zip(sequences, structures)
        ]
    elif backend == "threading":
        # Threading backend
        from concurrent.futures import ThreadPoolExecutor

        with ThreadPoolExecutor(max_workers=n_jobs) as executor:
            return list(
                executor.map(
                    lambda args: get_connectivity_list(args[0], args[1], format=format, **kwargs),
                    zip(sequences, structures),
                )
            )
    elif backend == "multiprocessing":
        # Multiprocessing backend
        with mp.Pool(processes=n_jobs) as pool:
            return pool.starmap(
                lambda seq, struct: get_connectivity_list(seq, struct, format=format, **kwargs),
                zip(sequences, structures),
            )
    else:
        raise ValueError(
            f"Unknown backend: {backend}. Must be 'multiprocessing', 'threading', or 'sequential'."
        )


def batch_apply(
    structs: List[SecStruct],
    func: Callable[[SecStruct], Any],
    n_jobs: Optional[int] = None,
    backend: str = "multiprocessing",
) -> List[Any]:
    """Apply a function to multiple structures in parallel.

    Args:
        structs: List of SecStruct objects.
        func: Function to apply to each structure.
        n_jobs: Number of parallel jobs. If None, uses CPU count.
        backend: Backend to use ('multiprocessing', 'threading', or 'sequential').

    Returns:
        List: List of function results.
    """
    if n_jobs is None:
        n_jobs = mp.cpu_count()

    if backend == "sequential" or n_jobs == 1:
        # Sequential processing
        return [func(s) for s in structs]
    elif backend == "threading":
        # Threading backend
        from concurrent.futures import ThreadPoolExecutor

        with ThreadPoolExecutor(max_workers=n_jobs) as executor:
            return list(executor.map(func, structs))
    elif backend == "multiprocessing":
        # Multiprocessing backend
        # Note: Functions must be picklable for multiprocessing
        with mp.Pool(processes=n_jobs) as pool:
            return pool.map(func, structs)
    else:
        raise ValueError(
            f"Unknown backend: {backend}. Must be 'multiprocessing', 'threading', or 'sequential'."
        )
