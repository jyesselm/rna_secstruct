# RNA Secondary Structure - Example Notebooks

This directory contains Jupyter notebooks demonstrating the functionality of the `rna_secstruct` library.

## Notebooks

1. **01_basic_usage.ipynb** - Basic usage, lazy loading, and fundamental operations
2. **02_connectivity.ipynb** - Connectivity lists, base pairs, and extended bracket notation
3. **03_structure_manipulation.ipynb** - Structure manipulation methods (insert, replace, remove, etc.)
4. **04_search_and_analysis.ipynb** - Search capabilities, pattern matching, and statistics
5. **05_json_serialization.ipynb** - JSON serialization and file I/O
6. **06_pandas_integration.ipynb** - Pandas integration with DataFrame/Series accessors
7. **07_parallel_processing.ipynb** - Parallel processing for batch operations

## Running the Notebooks

1. Install the required dependencies:
   ```bash
   pip install rna_secstruct pandas jupyter
   ```

2. Start Jupyter:
   ```bash
   jupyter notebook
   ```

3. Open and run the notebooks in order, or jump to specific topics of interest.

## Requirements

- Python 3.7+
- rna_secstruct (this package)
- pandas (for notebooks 06 and 07)
- jupyter (to run notebooks)

## Notes

- Notebooks are designed to be run independently, but following the order provides a logical progression
- Some notebooks require pandas (notebooks 06 and 07) - they will handle missing dependencies gracefully
- All notebooks include error handling examples and best practices

