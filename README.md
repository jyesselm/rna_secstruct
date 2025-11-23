# RNA Secondary Structure

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/downloads/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![License: Non-Commercial](https://img.shields.io/badge/license-Non--Commercial-yellow.svg)](LICENSE)

A modern Python package for parsing, analyzing, and manipulating RNA secondary structures. Designed with a clean API, lazy loading for performance, and comprehensive motif analysis capabilities.

## Features

✨ **Modern & Easy to Use** - Clean, intuitive API inspired by best practices  
🚀 **Performance Optimized** - Lazy loading for fast parsing of large structures  
🧬 **Comprehensive Analysis** - Extract, search, and manipulate structural motifs  
🔧 **Flexible Parsing** - Supports multiple bracket types, pseudoknots, and alternative formats  
📊 **Pandas Integration** - Seamless integration with pandas DataFrames  
⚡ **Parallel Processing** - Batch processing support for large datasets  
🛡️ **Robust Error Handling** - Graceful handling of malformed structures with warnings  

## Installation

Install from GitHub:

```bash
python -m pip install git+https://github.com/jyesselm/rna_secstruct
```

Install with optional dependencies:

```bash
# With pandas support
pip install git+https://github.com/jyesselm/rna_secstruct#egg=rna_secstruct[pandas]

# With parallel processing
pip install git+https://github.com/jyesselm/rna_secstruct#egg=rna_secstruct[parallel]

# With all optional dependencies
pip install git+https://github.com/jyesselm/rna_secstruct#egg=rna_secstruct[all]
```

## Quick Start

```python
from rna_secstruct import SecStruct

# Create a structure from sequence and dot-bracket notation
struct = SecStruct("GGGAAACCC", "(((...)))")

# Access basic properties
print(f"Sequence: {struct.sequence}")      # GGGAAACCC
print(f"Structure: {struct.structure}")    # (((...)))
print(f"Length: {len(struct)}")            # 9

# Access motifs (lazy loading - parsing happens here)
for motif_id, motif in struct.motifs.items():
    print(f"{motif_id}: {motif.m_type} - {motif.sequence}")
# 0: HELIX - GGG&CCC
# 1: HAIRPIN - GAAAC
```

## Examples

### Basic Usage

#### Creating Structures

```python
from rna_secstruct import SecStruct

# Simple hairpin
hairpin = SecStruct("GGGAAACCC", "(((...)))")

# Multi-strand structure (use & to separate strands)
multistrand = SecStruct(
    "GGGAAACCC&UUUAAA", 
    "(((...)))&(((...)))"
)

# Structure with junction
junction = SecStruct(
    "GGAAACGAAACGAAACC", 
    "((...)(...)(...))"
)
```

#### Accessing Motifs

```python
struct = SecStruct("GGGAAACCC", "(((...)))")

# Motifs are stored as a dictionary
print(struct.motifs)
# {0: HELIX,GGG&CCC,(((&))), 1: HAIRPIN,GAAAC,(...)}

# Access by ID
helix = struct[0]
hairpin = struct[1]

# Iterate over motifs
for motif in struct:
    print(f"{motif.m_type}: {motif.sequence}")
# HELIX: GGG&CCC
# HAIRPIN: GAAAC

# Get all motifs of a specific type
helices = struct.get_helices()
hairpins = struct.get_hairpins()
junctions = struct.get_junctions()
single_strands = struct.get_single_strands()
```

### Working with Motif Objects

```python
struct = SecStruct("GGGAAACCC", "(((...)))")
motif = struct[0]  # Get helix motif

# Basic properties
print(f"ID: {motif.m_id}")              # 0
print(f"Type: {motif.m_type}")          # HELIX
print(f"Sequence: {motif.sequence}")    # GGG&CCC
print(f"Structure: {motif.structure}")  # (((&)))

# Position information
print(f"Strands: {motif.strands}")      # [[0, 1, 2], [6, 7, 8]]
print(f"Positions: {motif.positions}")  # [0, 1, 2, 6, 7, 8]
print(f"Start: {motif.start_pos}")      # 0
print(f"End: {motif.end_pos}")          # 8

# Hierarchy
print(f"Has parent: {motif.has_parent()}")   # False
print(f"Has children: {motif.has_children()}")  # True
print(f"Children: {motif.children}")     # [HAIRPIN,GAAAC,(...)]

# Type checking
print(motif.is_helix())        # True
print(motif.is_hairpin())      # False
print(motif.is_junction())     # False
print(motif.is_single_strand())  # False

# Recursive operations (include all children)
seq, struct = motif.recursive_sequence(), motif.recursive_structure()
print(f"Recursive: {seq} {struct}")  # GGGAAACCC (((...)))
```

### Searching for Motifs

```python
from rna_secstruct import SecStruct, MotifSearchParams

struct = SecStruct("GGGACCUUCGGGACCC", "(((.((....)).)))")

# Search by sequence
results = struct.get_motifs(MotifSearchParams(sequence="GAC&GAC"))
print(results)
# [JUNCTION,GAC&GAC,(.(&).)]

# Search by structure pattern
results = struct.get_motifs(MotifSearchParams(structure="(....)"))
print(results)
# [HAIRPIN,GGAAAC,(....)]

# Search by motif type
helices = struct.get_motifs(MotifSearchParams(m_type="HELIX"))

# Search with position constraints (exclude 5' and 3' ends)
results = struct.get_motifs(
    MotifSearchParams(
        m_type="JUNCTION",
        min_pos=10,  # Start after position 10
        max_pos=50   # End before position 50
    )
)

# Search by token (motif identifier)
helix4 = struct.get_motifs_by_token("Helix4")  # Any helix of length 4
junction2 = struct.get_motifs_by_token("Junction2_5|0")  # 2-way junction with specific loop sizes
```

### Structure Manipulation

#### Changing Motifs

```python
struct = SecStruct("GGGAAACCC", "(((...)))")

# Change helix sequence
struct.change_motif(0, "AGG&CCU", "(((&)))")
print(struct.sequence)  # AGGAAACCU

# Change hairpin to hexaloop
struct.change_motif(1, "CUUUUUUG", "(......)")
print(struct.sequence)  # AGCUUUUUUGCU

# Replace with complex structure (auto-reparsing)
struct = SecStruct("GGGAAACCC", "(((...)))")
print("Before:")
print(struct.to_str())

struct.change_motif(1, "GGGACCUUCGGGACCC", "(((.((....)).)))")
print("\nAfter:")
print(struct.to_str())
# ID: 0, Helix5 GGGGG&CCCCC (((((&)))))
#    ID: 1, Junction2_1|1 GAC&GAC (.(&).)
#       ID: 2, Helix2 CC&GG ((&))
#          ID: 3, Hairpin4 CUUCGG (....)
```

#### Getting Substructures

```python
struct = SecStruct("GGGACCUUCGGGACCC", "(((.((....)).)))")

# Get a copy (important before making changes)
struct_copy = struct.get_copy()

# Get substructure starting from a motif
sub_struct = struct.get_sub_structure(1)  # From motif 1 and all its children
print(sub_struct.sequence)   # GACCUUCGGGAC
print(sub_struct.structure)  # (.((....)).)
```

### Connectivity Analysis

```python
from rna_secstruct import get_connectivity_list, ConnectivityList, STANDARD_BRACKET_TYPES

# Get connectivity list (pairmap)
struct = SecStruct("GGGAAACCC", "(((...)))")
conn = struct.connectivity
print(conn)  # [8, 7, 6, -1, -1, -1, 2, 1, 0]
# Index shows paired position, -1 means unpaired

# Check base pairs
cl = ConnectivityList("GGGAAACCC", "(((...)))")
print(cl.is_nucleotide_paired(0))    # True
print(cl.get_paired_nucleotide(0))   # 8
print(cl.get_basepair(0))            # GC

# Support for pseudoknots with multiple bracket types
pseudoknot = get_connectivity_list(
    "GGGAAACCC",
    "(([[))]]",
    bracket_types=STANDARD_BRACKET_TYPES  # Supports () [] {} <>
)
```

### Pandas Integration

```python
import pandas as pd
from rna_secstruct import SecStruct

# Create a DataFrame with sequences and structures
df = pd.DataFrame({
    'sequence': ['GGGAAACCC', 'GGAAACGAAAC', 'GGGACCUUCGGGACCC'],
    'structure': ['(((...)))', '((...)(...))', '(((.((....)).)))']
})

# Convert to SecStruct objects
df['secstruct'] = df.apply(
    lambda row: SecStruct(row['sequence'], row['structure']), 
    axis=1
)

# Access motifs directly
df['num_helices'] = df['secstruct'].apply(lambda s: len(s.get_helices()))
df['num_hairpins'] = df['secstruct'].apply(lambda s: len(s.get_hairpins()))

# Or use the accessor (if registered)
df['secstruct'].secstruct.get_helices()  # Returns list of lists
```

### Parallel Processing

```python
from rna_secstruct import batch_parse
import pandas as pd

# Large dataset
sequences = ["GGGAAACCC"] * 1000
structures = ["(((...)))"] * 1000

# Process in parallel
results = batch_parse(sequences, structures, n_jobs=4)

# Or use pandas extension
df = pd.DataFrame({
    'sequence': sequences,
    'structure': structures
})
results = df.apply(
    lambda row: SecStruct(row['sequence'], row['structure']),
    axis=1
)
```

### Working with Real-World Data

```python
from rna_secstruct import SecStruct, MotifSearchParams

# Large RNA structure
seq = (
    "GGAAGAUCGAGUAGAUCAAAGAGCCUAUGGCUGCCACCCGAGCCCUUGAACUACAGGGAACACUGGAAA"
    "CAGUACCCCCUGCAAGGGCGUUUGACGGUGGCAGCCUAAGGGCUCAAAGAAACAACAACAACAAC"
)
ss = (
    "....((((.....))))...((((((..((((((((((((((((((((.....(((((...((((....)"
    ")))...))))))))))))..)))..))))))))))...))))))...................."
)

struct = SecStruct(seq, ss)

# Find all junctions after position 50
junctions = struct.get_motifs(
    MotifSearchParams(m_type="JUNCTION", min_pos=50)
)

# Find all 5-nucleotide hairpins
hairpins_5 = struct.get_motifs(MotifSearchParams(structure="(....)"))

# Get motif statistics
print(f"Total motifs: {len(struct.motifs)}")
print(f"Helices: {len(struct.get_helices())}")
print(f"Hairpins: {len(struct.get_hairpins())}")
print(f"Junctions: {len(struct.get_junctions())}")
```

### Error Handling

The parser handles invalid inputs gracefully with warnings:

```python
import logging
from rna_secstruct import Parser

# Set up logging to see warnings
logging.basicConfig(level=logging.WARNING)

p = Parser()

# These will log warnings but still parse:
# - Invalid characters (replaced with 'N' or '.')
# - Length mismatches (truncated/padded)
# - Unbalanced parentheses (auto-balanced)
# - Invalid bracket types (normalized)

result = p.parse("GGGAAACCC", "(((...)))(")  # Unbalanced - will auto-fix
result = p.parse("GGGYAACCC", "(((...)))")   # Invalid 'Y' - replaced with 'N'
result = p.parse("GGGAAACCC", "((([...)))")  # Invalid bracket - normalized
```

### Advanced: Multi-Strand Structures

```python
from rna_secstruct import SecStruct

# Two separate RNA molecules
struct = SecStruct(
    "GGGAAACCC&UUUGGGAAA", 
    "(((...)))&(((...)))"
)

# Access strands separately
print(struct.sequence.count('&'))  # Number of strand separators

# Iterate over motifs (includes all strands)
for motif in struct:
    print(motif.sequence)  # May contain '&' for multi-strand motifs
```

### Advanced: Pseudoknot Support

```python
from rna_secstruct import SecStruct, STANDARD_BRACKET_TYPES

# Pseudoknot structure using different bracket types
pseudoknot = SecStruct("GGGAAACCC", "(([[))]]")

# The parser preserves bracket types for pseudoknot representation
# Use connectivity module for full pseudoknot analysis
from rna_secstruct import get_connectivity_list

conn = get_connectivity_list(
    "GGGAAACCC",
    "(([[))]]",
    bracket_types=STANDARD_BRACKET_TYPES  # () [] {} <>
)
```

## API Overview

### Main Classes

- **`SecStruct`** - Main class for RNA secondary structures
- **`Motif`** - Represents individual structural motifs
- **`MotifSearchParams`** - Parameters for motif searching
- **`ConnectivityList`** - Connectivity/pairmap representation

### Key Methods

#### SecStruct Methods
- `get_motifs(params)` - Search for motifs with constraints
- `get_motifs_by_token(token)` - Search by motif identifier
- `get_helices()`, `get_hairpins()`, `get_junctions()` - Get specific motif types
- `change_motif(id, sequence, structure)` - Modify a motif
- `get_sub_structure(id)` - Extract substructure
- `get_copy()` - Create a copy
- `to_str()` - Format structure representation

#### Motif Properties
- `m_id`, `m_type`, `sequence`, `structure`
- `strands`, `positions`, `start_pos`, `end_pos`
- `parent`, `children`
- `recursive_sequence()`, `recursive_structure()`

## Documentation

- **Jupyter Notebooks**: See `notebooks/` directory for detailed examples
- **API Documentation**: Check docstrings in source code
- **Examples**: All examples in this README are runnable

## Development

### Running Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=rna_secstruct --cov-report=html

# Run specific test file
pytest test/test_parser.py
```

### Code Quality

```bash
# Format code
black rna_secstruct/ test/

# Lint
ruff check rna_secstruct/ test/

# Type checking
mypy rna_secstruct/
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under a Non-Commercial License. See [LICENSE](LICENSE) file for details.

For commercial licensing inquiries, please contact: jyesselm@unl.edu

## Citation

If you use `rna_secstruct` in your research, please cite:

```bibtex
@software{rna_secstruct,
  author = {Yesselman, Joe},
  title = {rna_secstruct: A Python package for RNA secondary structure analysis},
  url = {https://github.com/jyesselm/rna_secstruct},
  version = {0.1.0},
  year = {2024}
}
```

## Links

- **GitHub**: https://github.com/jyesselm/rna_secstruct
- **Issues**: https://github.com/jyesselm/rna_secstruct/issues
- **Author**: Joe Yesselman (jyesselm@unl.edu)

---

**Note**: This package is designed for non-commercial use. For commercial applications, please contact the author for licensing options.
