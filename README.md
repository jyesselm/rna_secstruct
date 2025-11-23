# RNA Secondary Structure

[![PyPI version](https://badge.fury.io/py/rna-secstruct.svg)](https://badge.fury.io/py/rna-secstruct)
[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://github.com/jyesselm/rna_secstruct/actions/workflows/tests.yml/badge.svg)](https://github.com/jyesselm/rna_secstruct/actions/workflows/tests.yml)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![License: Non-Commercial](https://img.shields.io/badge/license-Non--Commercial-yellow.svg)](LICENSE)

A modern, comprehensive Python package for parsing, analyzing, and manipulating RNA secondary structures. Designed with a clean API, lazy loading for performance, comprehensive motif analysis, and extensive integration capabilities.

## Features

✨ **Modern & Easy to Use** - Clean, intuitive API inspired by best practices  
🚀 **Performance Optimized** - Lazy loading for fast parsing of large structures  
🧬 **Comprehensive Analysis** - Extract, search, and manipulate structural motifs  
🔧 **Flexible Parsing** - Supports multiple bracket types, pseudoknots, and alternative formats  
📊 **Pandas Integration** - Seamless integration with pandas DataFrames via accessors  
⚡ **Parallel Processing** - Batch processing support for large datasets  
🛡️ **Robust Error Handling** - Graceful handling of malformed structures with warnings  
🔍 **Type Safe** - Full type annotations with mypy support for better code quality  
📦 **JSON Serialization** - Built-in JSON support for data exchange  
🔬 **Advanced Search** - Pattern matching, wildcards, and complex motif queries  
🔄 **Structure Manipulation** - Immutable operations for safe structure modification  
📈 **Statistics & Analysis** - Comprehensive metrics and comparison tools  
✅ **Validation** - Built-in validation and normalization utilities  

## Installation

Install from PyPI:

```bash
pip install rna_secstruct
```

Install with optional dependencies:

```bash
# With pandas support
pip install rna_secstruct[pandas]

# With parallel processing
pip install rna_secstruct[parallel]

# With all optional dependencies
pip install rna_secstruct[all]
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

## Table of Contents

- [Basic Usage](#basic-usage)
- [Working with Motifs](#working-with-motifs)
- [Advanced Search](#advanced-search)
- [Structure Manipulation](#structure-manipulation)
- [Connectivity Analysis](#connectivity-analysis)
- [Statistics & Analysis](#statistics--analysis)
- [Validation & Normalization](#validation--normalization)
- [Comparison Operations](#comparison-operations)
- [JSON Serialization](#json-serialization)
- [Pandas Integration](#pandas-integration)
- [Parallel Processing](#parallel-processing)
- [Multi-Strand Structures](#multi-strand-structures)
- [Pseudoknot Support](#pseudoknot-support)
- [Error Handling](#error-handling)
- [API Reference](#api-reference)

## Basic Usage

### Creating Structures

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

# Complex structure
complex_struct = SecStruct(
    "GGGACCUUCGGGACCC",
    "(((.((....)).)))"
)
```

### Basic Properties

```python
struct = SecStruct("GGGAAACCC", "(((...)))")

# Sequence and structure
print(struct.sequence)    # GGGAAACCC
print(struct.structure)   # (((...)))
print(len(struct))        # 9

# String representation
print(repr(struct))       # GGGAAACCC, (((...)))

# Connectivity (pairmap)
print(struct.connectivity)  # [8, 7, 6, -1, -1, -1, 2, 1, 0]
# -1 means unpaired, numbers indicate paired position
```

### Slicing and Indexing

```python
struct = SecStruct("GGGAAACCC", "(((...)))")

# Slicing returns new SecStruct (immutable)
substruct = struct[3:6]  # Extract middle region
print(substruct.sequence)   # AAA
print(substruct.structure)  # ...

# Access motifs by ID (after parsing)
motif = struct[0]  # Get motif with ID 0
print(motif.m_type)  # HELIX
```

### Concatenation

```python
struct1 = SecStruct("GGG", "(((")
struct2 = SecStruct("AAA", "...")
combined = struct1 + struct2  # New SecStruct instance
print(combined.sequence)  # GGGAAA
```

## Working with Motifs

### Accessing Motifs

```python
struct = SecStruct("GGGAAACCC", "(((...)))")

# Motifs are stored as a dictionary (lazy-loaded)
print(struct.motifs)
# {0: HELIX,GGG&CCC,(((&))), 1: HAIRPIN,GAAAC,(...)}

# Access by ID
helix = struct[0]
hairpin = struct[1]

# Iterate over motifs
for motif in struct:
    print(f"{motif.m_type}: {motif.sequence}")

# Iterate with IDs
for motif_id, motif in struct.itermotifs():
    print(f"ID {motif_id}: {motif.m_type}")

# Get specific motif types
helices = struct.get_helices()
hairpins = struct.get_hairpins()
junctions = struct.get_junctions()
single_strands = struct.get_single_strands()

# Count motifs
print(struct.get_num_motifs())  # Total number of motifs
```

### Motif Properties

```python
struct = SecStruct("GGGAAACCC", "(((...)))")
motif = struct[0]  # Get helix motif

# Basic properties
print(f"ID: {motif.m_id}")              # 0
print(f"Type: {motif.m_type}")          # HELIX
print(f"Sequence: {motif.sequence}")     # GGG&CCC
print(f"Structure: {motif.structure}")   # (((&)))
print(f"Token: {motif.token}")          # Helix3

# Position information
print(f"Strands: {motif.strands}")      # [[0, 1, 2], [6, 7, 8]]
print(f"Positions: {motif.positions}")    # [0, 1, 2, 6, 7, 8]
print(f"Start: {motif.start_pos}")      # 0
print(f"End: {motif.end_pos}")          # 8

# Hierarchy
print(f"Has parent: {motif.has_parent()}")      # False
print(f"Has children: {motif.has_children()}")   # True
print(f"Children: {motif.children}")              # [HAIRPIN,GAAAC,(...)]
print(f"Parent: {motif.parent}")                  # None

# Type checking
print(motif.is_helix())         # True
print(motif.is_hairpin())       # False
print(motif.is_junction())      # False
print(motif.is_single_strand()) # False

# Position checking
print(motif.contains(5))        # True (position 5 is in this motif)

# Recursive operations (include all children)
seq = motif.recursive_sequence()
struct = motif.recursive_structure()
print(f"Recursive: {seq} {struct}")  # GGGAAACCC (((...)))

# String representation
print(motif.to_str())
# ID: 0, Helix3 GGG&CCC (((&)))
#    ID: 1, Hairpin5 GAAAC (...)
```


## Advanced Search

### Basic Motif Search

```python
from rna_secstruct import SecStruct, MotifSearchParams

struct = SecStruct("GGGACCUUCGGGACCC", "(((.((....)).)))")

# Search by sequence
results = struct.get_motifs(MotifSearchParams(sequence="GAC&GAC"))
print(results)  # [JUNCTION,GAC&GAC,(.(&).)]

# Search by structure pattern
results = struct.get_motifs(MotifSearchParams(structure="(....)"))
print(results)  # [HAIRPIN,GGAAAC,(....)]

# Search by motif type
helices = struct.get_motifs(MotifSearchParams(m_type="HELIX"))
hairpins = struct.get_motifs(MotifSearchParams(m_type="HAIRPIN"))
junctions = struct.get_motifs(MotifSearchParams(m_type="JUNCTION"))
```

### Advanced Search Parameters

```python
struct = SecStruct("GGGACCUUCGGGACCC", "(((.((....)).)))")

# Position constraints
results = struct.get_motifs(
    MotifSearchParams(
        m_type="JUNCTION",
        min_pos=10,  # Start after position 10
        max_pos=50    # End before position 50
    )
)

# Length constraints
results = struct.get_motifs(
    MotifSearchParams(
        m_type="HAIRPIN",
        min_length=4,   # At least 4 nucleotides
        max_length=10   # At most 10 nucleotides
    )
)

# ID constraints
results = struct.get_motifs(
    MotifSearchParams(
        min_id=1,      # Motif ID >= 1
        max_id=5       # Motif ID <= 5
    )
)

# Children constraints
results = struct.get_motifs(
    MotifSearchParams(
        m_type="HELIX",
        has_children=True  # Only helices with children
    )
)

# Combined search
results = struct.get_motifs(
    MotifSearchParams(
        m_type="HAIRPIN",
        min_length=4,
        max_length=8,
        min_pos=5,
        max_pos=20
    )
)
```

### Token-Based Search

```python
struct = SecStruct("GGGAAACCC", "(((...)))")

# Search by token (motif identifier)
helix4 = struct.get_motifs_by_token("Helix4")  # Any helix of length 4
junction2 = struct.get_motifs_by_token("Junction2_5|0")  # 2-way junction

# Token format examples:
# - "Helix3" - Helix with 3 base pairs
# - "Hairpin5" - Hairpin with 5 nucleotides
# - "Junction2_3|4" - 2-way junction with loop sizes 3 and 4
```

### Strand Length Search

```python
struct = SecStruct("GGGACCUUCGGGACCC", "(((.((....)).)))")

# Find motifs by strand lengths
# For a hairpin: [5] means 5 nucleotides
hairpins_5 = struct.get_motifs_by_strand_lengths([5])

# For a junction: [3, 3] means two strands of length 3
junctions_3_3 = struct.get_motifs_by_strand_lengths([3, 3])

# For a helix: [3, 3] means two strands of length 3
helices_3_3 = struct.get_motifs_by_strand_lengths([3, 3])
```

### Topology-Based Junction Search

```python
struct = SecStruct("GGGACCUUCGGGACCC", "(((.((....)).)))")

# Find two-way junctions by topology
# x_pos and y_pos are loop sizes (excluding flanking base pairs)
junctions = struct.get_twoway_junctions_by_topology(x_pos=1, y_pos=1)
# This finds junctions with loop sizes 1 and 1
```

### Sequence and Structure Pattern Matching

```python
struct = SecStruct("GGGAAACCC", "(((...)))")

# Find sequence patterns (with wildcards)
matches = struct.find_sequence("GAA", allow_wildcards=False)
# Returns: [(3, 6)]  # (start, end) positions

# Find with wildcards
matches = struct.find_sequence("GNN", allow_wildcards=True)
# N matches any nucleotide
# R = A/G, Y = U/C, M = A/C, K = U/G, S = G/C, W = A/U
# B = not A, D = not C, H = not G, V = not U

# Find structure patterns
matches = struct.find_structure("(((")
# Returns: [(0, 3)]  # (start, end) positions

# Find complete substructures
sub = SecStruct("AAA", "...")
matches = struct.find(sub)
# Returns: [(3, 6)]  # (start, end) positions
```

## Structure Manipulation

Most manipulation operations return **new** `SecStruct` instances (immutable pattern). The original structure is never modified. However, `change_motif()` modifies the structure in place.

### Changing Motifs

```python
struct = SecStruct("GGGAAACCC", "(((...)))")

# Change helix sequence
# Note: change_motif() modifies the SecStruct in place (not immutable)
struct.change_motif(0, "AGG&CCU", "(((&)))")
print(struct.sequence)  # AGGAAACCU

# Change hairpin to hexaloop
struct.change_motif(1, "CUUUUUUG", "(......)")
print(struct.sequence)  # AGGCUUUUUUGCCU

# Replace with complex structure (auto-reparsing)
struct = SecStruct("GGGAAACCC", "(((...)))")
struct.change_motif(1, "GGGACCUUCGGGACCC", "(((.((....)).)))")
print(struct.to_str())
# ID: 0, Helix5 GGGGG&CCCCC (((((&)))))
#    ID: 1, Junction2_1|1 GAC&GAC (.(&).)
#       ID: 2, Helix2 CC&GG ((&))
#          ID: 3, Hairpin4 CUUCGG (....)
```

### Getting Substructures

```python
struct = SecStruct("GGGACCUUCGGGACCC", "(((.((....)).)))")

# Get a copy (important before making changes)
struct_copy = struct.get_copy()

# Get substructure starting from a motif
sub_struct = struct.get_sub_structure(1)  # From motif 1 and all its children
print(sub_struct.sequence)   # GACCUUCGGGAC
print(sub_struct.structure)  # (.((....)).)
```

### Immutable Container Operations

```python
struct1 = SecStruct("GGG", "(((")
struct2 = SecStruct("AAA", "...")
struct3 = SecStruct("CCC", ")))")

# Split strands (if multi-strand)
strands = struct.split_strands()  # Returns list of SecStruct objects

# Insert at position
new_struct = struct1.insert(3, struct2)
# Returns: SecStruct("GGGAAA", "(((...)")

# Join structures (with & separator)
joined = struct1.join(struct2)
# Returns: SecStruct("GGG&AAA", "(((&...)")

# Replace at position
replaced = struct1.replace(struct2, 0)
# Returns: SecStruct("AAA", "...")

# Remove region
removed = struct1.remove(1, 2)
# Returns: SecStruct("GG", "((")

# Subtract substructure
subtracted = struct1.subtract(struct2)
# Removes struct2 from struct1 if found
```

## Connectivity Analysis

### Basic Connectivity

```python
from rna_secstruct import SecStruct, get_connectivity_list, ConnectivityList

# Get connectivity list (pairmap)
struct = SecStruct("GGGAAACCC", "(((...)))")
conn = struct.connectivity
print(conn)  # [8, 7, 6, -1, -1, -1, 2, 1, 0]
# Index shows paired position, -1 means unpaired

# Using ConnectivityList class
cl = ConnectivityList("GGGAAACCC", "(((...)))")
print(cl.connections)  # [8, 7, 6, -1, -1, -1, 2, 1, 0]
print(cl.sequence)    # GGGAAACCC
print(cl.structure)   # (((...)))
```

### Base Pair Operations

```python
cl = ConnectivityList("GGGAAACCC", "(((...)))")

# Check if nucleotide is paired
print(cl.is_nucleotide_paired(0))    # True
print(cl.is_nucleotide_paired(3))    # False

# Get paired nucleotide
print(cl.get_paired_nucleotide(0))   # 8
print(cl.get_paired_nucleotide(8))   # 0

# Get base pair
print(cl.get_basepair(0))            # GC
print(cl.get_basepair(3))            # . (unpaired)

# Get pair type (bracket character, letter, or number)
print(cl.get_pair_type(0))           # '('
```

### Using SecStruct Connectivity Methods

```python
struct = SecStruct("GGGAAACCC", "(((...)))")

# Check if position is paired
print(struct.is_paired(0))           # True

# Get base pair tuple
bp = struct.get_basepair(0)
print(bp)                            # (0, 8) or None if unpaired

# Count base pairs
print(struct.get_num_basepairs())     # 3

# Count unpaired
print(struct.get_num_unpaired())     # 3
```

### Pseudoknot Support

```python
from rna_secstruct import get_connectivity_list, STANDARD_BRACKET_TYPES

# Pseudoknot structure using different bracket types
pseudoknot = get_connectivity_list(
    "GGGAAACCC",
    "(([[))]]",
    bracket_types=STANDARD_BRACKET_TYPES  # Supports () [] {} <>
)

print(pseudoknot.is_nucleotide_paired(0))  # True
print(pseudoknot.get_pair_type(0))         # '('
print(pseudoknot.get_pair_type(3))         # '['

# Detect pseudoknots
from rna_secstruct.connectivity import has_pseudoknot
conn = pseudoknot.connections
has_pk = has_pseudoknot(conn, STANDARD_BRACKET_TYPES)
print(has_pk)  # True
```

### Alternative Formats

```python
from rna_secstruct import get_connectivity_list

# Letter-based format
cl = get_connectivity_list(
    "GGGAAACCC",
    "a b c c b a",
    format="letter"
)
print(cl.get_pair_type(0))  # 'a'

# Number-based format
cl = get_connectivity_list(
    "GGGAAACCC",
    "1 2 3 3 2 1",
    format="number"
)
print(cl.get_pair_type(0))  # '1'

# Auto-detect format
cl = get_connectivity_list(
    "GGGAAACCC",
    "(((...)))"
    # format=None auto-detects
)
```

### Circular Structure Detection

```python
from rna_secstruct.connectivity import is_circular

conn = [1, 2, 0]  # Circular: 0->1, 1->2, 2->0
is_circ = is_circular(0, conn)
print(is_circ)  # True
```

## Statistics & Analysis

### Basic Statistics

```python
struct = SecStruct("GGGAAACCC", "(((...)))")

# Base pair statistics
print(struct.get_num_basepairs())    # 3
print(struct.get_num_unpaired())    # 3

# GC content
print(struct.get_gc_content())       # 0.666... (2/3)

# Helix information
helix_lengths = struct.get_helix_lengths()
print(helix_lengths)                # [3] (lengths of all helices)

# Motif counts
print(struct.get_num_motifs())       # 2
print(len(struct.get_helices()))     # 1
print(len(struct.get_hairpins()))    # 1
print(len(struct.get_junctions()))   # 0
```

### Structure Comparison

```python
struct1 = SecStruct("GGGAAACCC", "(((...)))")
struct2 = SecStruct("GGGAAACCC", "(((...)))")
struct3 = SecStruct("AAA", "...")

# Equality
print(struct1 == struct2)            # True
print(struct1 == struct3)           # False

# Structural similarity (structure string comparison)
similarity = struct1.structural_similarity(struct2)
print(similarity)                    # 1.0 (identical structures)

# Sequence identity
identity = struct1.sequence_identity(struct2)
print(identity)                      # 1.0 (identical sequences)
```

## Validation & Normalization

### Validation

```python
struct = SecStruct("GGGAAACCC", "(((...)))")

# Validate structure (raises ValueError if invalid)
try:
    struct.validate()
    print("Structure is valid!")
except ValueError as e:
    print(f"Invalid structure: {e}")

# Check validity (non-raising)
if struct.is_valid():
    print("Structure is valid")
else:
    print("Structure is invalid")
```

### Normalization

```python
struct = SecStruct("gggaaaccc", "(((...)))")

# Normalize (uppercase, T->U conversion)
normalized = struct.normalize()
print(normalized.sequence)  # GGGAAACCC (uppercase)
# T nucleotides are converted to U
```

## Comparison Operations

```python
struct1 = SecStruct("GGGAAACCC", "(((...)))")
struct2 = SecStruct("GGGAAACCC", "(((...)))")
struct3 = SecStruct("AAA", "...")

# Equality
print(struct1 == struct2)  # True
print(struct1 == struct3)  # False

# Structural similarity (0.0 to 1.0)
similarity = struct1.structural_similarity(struct2)
print(similarity)  # 1.0

# Sequence identity (0.0 to 1.0)
identity = struct1.sequence_identity(struct2)
print(identity)  # 1.0
```

## JSON Serialization

### Basic JSON Operations

```python
from rna_secstruct import SecStruct

struct = SecStruct("GGGAAACCC", "(((...)))")

# Convert to dictionary
data = struct.to_dict()
print(data)
# {'sequence': 'GGGAAACCC', 'structure': '(((...)))', 'motifs': [...]}

# Convert to JSON string
json_str = struct.to_json(indent=2)
print(json_str)

# Create from dictionary
struct2 = SecStruct.from_dict(data)

# Create from JSON string
struct3 = SecStruct.from_json(json_str)
```

### File Operations

```python
struct = SecStruct("GGGAAACCC", "(((...)))")

# Save to file
struct.to_json_file("structure.json", indent=2)

# Load from file
loaded = SecStruct.from_json_file("structure.json")
```

### Custom JSON Encoder

```python
from rna_secstruct.json_encoder import SecStructJSONEncoder, dumps, loads
import json

struct = SecStruct("GGGAAACCC", "(((...)))")

# Use custom encoder
json_str = json.dumps(struct, cls=SecStructJSONEncoder, indent=2)

# Or use convenience function
json_str = dumps(struct, indent=2)

# Load back
data = loads(json_str)
struct2 = SecStruct.from_dict(data)
```

### Comma-Delimited Format

```python
struct = SecStruct("GGGAAACCC", "(((...)))")

# CSV representation
csv = struct.to_comma_delimited()
print(csv)  # GGGAAACCC,(((...)))
```

## Pandas Integration

### Basic Usage

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
```

### DataFrame Accessor

```python
import pandas as pd
from rna_secstruct import SecStruct

df = pd.DataFrame({
    'sequence': ['GGGAAACCC', 'GGAAACGAAAC'],
    'structure': ['(((...)))', '((...)(...))']
})

# Create SecStruct column using accessor
df = df.rna.add_secstruct('sequence', 'structure', column='secstruct')

# Add statistics columns
df = df.rna.add_statistics('secstruct')
# Adds: secstruct_num_bp, secstruct_num_unpaired, 
#       secstruct_gc_content, secstruct_length
```

### Series Accessor

```python
import pandas as pd
from rna_secstruct import SecStruct

# Series of SecStruct objects
series = pd.Series([
    SecStruct("GGGAAACCC", "(((...)))"),
    SecStruct("GGAAACGAAAC", "((...)(...))")
])

# Get statistics
num_bp = series.rna.num_basepairs()
num_motifs = series.rna.num_motifs()
gc_content = series.rna.gc_content()
helix_lengths = series.rna.helix_lengths()
has_pk = series.rna.has_pseudoknot()

# JSON operations
json_str = series.rna.to_json(indent=2)
series2 = series.rna.from_json(json_str)
```

## Parallel Processing

### Batch Parsing

```python
from rna_secstruct import batch_parse

# Large dataset
sequences = ["GGGAAACCC"] * 1000
structures = ["(((...)))"] * 1000

# Process in parallel
results = batch_parse(
    sequences, 
    structures, 
    n_jobs=4,  # Number of parallel jobs
    backend="multiprocessing"  # or "threading" or "sequential"
)

print(len(results))  # 1000
```

### Batch Connectivity

```python
from rna_secstruct.parallel import batch_connectivity

sequences = ["GGGAAACCC"] * 1000
structures = ["(((...)))"] * 1000

# Generate connectivity lists in parallel
conn_lists = batch_connectivity(
    sequences,
    structures,
    n_jobs=4,
    backend="multiprocessing"
)
```

### Batch Apply

```python
from rna_secstruct.parallel import batch_apply
from rna_secstruct import SecStruct

# List of structures
structs = [SecStruct("GGGAAACCC", "(((...)))")] * 1000

# Apply function in parallel
def count_motifs(s):
    return s.get_num_motifs()

results = batch_apply(
    structs,
    count_motifs,
    n_jobs=4,
    backend="multiprocessing"
)
```

### Backend Options

```python
# Multiprocessing (default, good for CPU-bound tasks)
results = batch_parse(seqs, structs, backend="multiprocessing", n_jobs=4)

# Threading (good for I/O-bound tasks)
results = batch_parse(seqs, structs, backend="threading", n_jobs=4)

# Sequential (no parallelization)
results = batch_parse(seqs, structs, backend="sequential")
```

## Multi-Strand Structures

### Creating Multi-Strand Structures

```python
from rna_secstruct import SecStruct

# Two separate RNA molecules
struct = SecStruct(
    "GGGAAACCC&UUUGGGAAA", 
    "(((...)))&(((...)))"
)

# Access strands separately
print(struct.sequence.count('&'))  # Number of strand separators

# Split into individual strands
strands = struct.split_strands()
for i, strand in enumerate(strands):
    print(f"Strand {i}: {strand.sequence} {strand.structure}")

# Iterate over motifs (includes all strands)
for motif in struct:
    print(motif.sequence)  # May contain '&' for multi-strand motifs
```

### Multi-Strand Operations

```python
struct1 = SecStruct("GGG", "(((")
struct2 = SecStruct("AAA", "...")

# Join with & separator
joined = struct1.join(struct2)
print(joined.sequence)  # GGG&AAA
```

## Pseudoknot Support

### Basic Pseudoknots

```python
from rna_secstruct import SecStruct, get_connectivity_list, STANDARD_BRACKET_TYPES

# Pseudoknot structure using different bracket types
pseudoknot = SecStruct("GGGAAACCC", "(([[))]]")

# Use connectivity module for full pseudoknot analysis
conn = get_connectivity_list(
    "GGGAAACCC",
    "(([[))]]",
    bracket_types=STANDARD_BRACKET_TYPES  # Supports () [] {} <>
)

print(conn.is_nucleotide_paired(0))  # True
print(conn.get_pair_type(0))         # '('
print(conn.get_pair_type(3))         # '['
```

### Detecting Pseudoknots

```python
from rna_secstruct.connectivity import has_pseudoknot, STANDARD_BRACKET_TYPES

# Simple structure (no pseudoknot)
conn1 = [8, 7, 6, -1, -1, -1, 2, 1, 0]
print(has_pseudoknot(conn1, STANDARD_BRACKET_TYPES))  # False

# Pseudoknot structure
conn2 = get_connectivity_list("GGGAAACCC", "(([[))]]", 
                               bracket_types=STANDARD_BRACKET_TYPES)
print(has_pseudoknot(conn2.connections, STANDARD_BRACKET_TYPES))  # True
```

## Error Handling

The parser handles invalid inputs gracefully with warnings:

```python
import logging
from rna_secstruct import Parser

# Set up logging to see warnings
logging.basicConfig(level=logging.WARNING)

p = Parser()

# These will log warnings but still parse:
# - Invalid characters (replaced with 'N' or '.')
result = p.parse("GGGYAACCC", "(((...)))")   # Invalid 'Y' - replaced with 'N'

# - Length mismatches (truncated/padded)
result = p.parse("GGGAAACCC", "(((...)))(")   # Unbalanced - will auto-fix

# - Unbalanced parentheses (auto-balanced)
result = p.parse("GGGAAACCC", "((([...)))")   # Invalid bracket - normalized

# - Invalid bracket types (normalized)
result = p.parse("GGGAAACCC", "(((...)))")   # Valid structure
```

### Validation Errors

```python
from rna_secstruct import SecStruct

# These will raise ValueError:
try:
    # Length mismatch
    struct = SecStruct("GGG", "(((")  # OK
    struct = SecStruct("GGG", "(((")  # OK
except ValueError as e:
    print(f"Error: {e}")

# Invalid structure (if validation enabled)
try:
    struct = SecStruct("GGGAAACCC", "(((...)))")
    struct.validate()  # Raises if invalid
except ValueError as e:
    print(f"Validation error: {e}")
```

## API Reference

### Main Classes

#### `SecStruct`

Main class for RNA secondary structures.

**Key Methods:**
- `get_motifs(params)` - Search for motifs with constraints
- `get_motifs_by_token(token)` - Search by motif identifier
- `get_motifs_by_strand_lengths(lengths)` - Search by strand lengths
- `get_twoway_junctions_by_topology(x, y)` - Find junctions by topology
- `get_helices()`, `get_hairpins()`, `get_junctions()`, `get_single_strands()` - Get specific motif types
- `change_motif(id, sequence, structure)` - Modify a motif
- `get_sub_structure(id)` - Extract substructure
- `get_copy()` - Create a copy
- `to_str()` - Format structure representation
- `split_strands()` - Split multi-strand structure
- `insert(pos, other)` - Insert structure at position
- `join(other)` - Join structures with &
- `replace(other, pos)` - Replace at position
- `remove(start, end)` - Remove region
- `subtract(other)` - Remove substructure
- `find(sub)` - Find substructure positions
- `find_sequence(pattern)` - Find sequence pattern
- `find_structure(pattern)` - Find structure pattern
- `is_paired(index)` - Check if position is paired
- `get_basepair(index)` - Get base pair tuple
- `get_num_basepairs()` - Count base pairs
- `get_num_unpaired()` - Count unpaired nucleotides
- `get_gc_content()` - Calculate GC content
- `get_helix_lengths()` - Get helix lengths
- `validate()` - Validate structure
- `is_valid()` - Check validity
- `normalize()` - Normalize sequence
- `structural_similarity(other)` - Compare structures
- `sequence_identity(other)` - Compare sequences
- `to_dict()` - Convert to dictionary
- `to_json()` - Serialize to JSON
- `from_dict(data)` - Create from dictionary
- `from_json(json_str)` - Deserialize from JSON
- `to_json_file(filepath)` - Save to file
- `from_json_file(filepath)` - Load from file

**Properties:**
- `sequence` - RNA sequence
- `structure` - Secondary structure
- `motifs` - Dictionary of motifs (lazy-loaded)
- `connectivity` - Connectivity list (pairmap)

#### `Motif`

Represents individual structural motifs.

**Properties:**
- `m_id` - Motif ID
- `m_type` - Motif type (HELIX, HAIRPIN, JUNCTION, SINGLESTRAND)
- `sequence` - Motif sequence
- `structure` - Motif structure
- `strands` - List of strand indices
- `positions` - All positions in motif
- `start_pos` - Start position
- `end_pos` - End position
- `parent` - Parent motif
- `children` - List of child motifs
- `token` - Motif identifier token

**Methods:**
- `contains(position)` - Check if position is in motif
- `has_parent()` - Check if has parent
- `has_children()` - Check if has children
- `is_helix()`, `is_hairpin()`, `is_junction()`, `is_single_strand()` - Type checks
- `num_strands()` - Number of strands
- `recursive_sequence()` - Sequence including children
- `recursive_structure()` - Structure including children
- `to_str(depth)` - String representation
- `to_dict()` - Convert to dictionary

#### `MotifSearchParams`

Parameters for RNA motif search.

**Attributes:**
- `sequence` - Exact sequence to match
- `structure` - Exact structure to match
- `m_type` - Motif type to match
- `min_pos`, `max_pos` - Position constraints
- `min_id`, `max_id` - ID constraints
- `token` - Token to match
- `min_length`, `max_length` - Length constraints
- `strand_lengths` - List of strand lengths to match
- `has_children` - Whether motif must have children

#### `ConnectivityList`

Connectivity/pairmap representation.

**Methods:**
- `is_nucleotide_paired(index)` - Check if paired
- `get_paired_nucleotide(index)` - Get paired position
- `get_basepair(index)` - Get base pair string
- `get_pair_type(index)` - Get pair type (bracket/letter/number)

**Properties:**
- `connections` - Connectivity list
- `sequence` - RNA sequence
- `structure` - Secondary structure
- `pair_types` - Dictionary of pair types

### Utility Functions

- `get_connectivity_list(sequence, structure, format, bracket_types)` - Create ConnectivityList
- `connectivity_list(structure, bracket_types)` - Get simple connectivity list
- `has_pseudoknot(connectivity_lists, bracket_types)` - Detect pseudoknots
- `is_circular(start, connections)` - Detect circular structures
- `batch_parse(sequences, structures, n_jobs, backend)` - Parallel parsing
- `batch_connectivity(sequences, structures, format, n_jobs, backend)` - Parallel connectivity
- `batch_apply(structs, func, n_jobs, backend)` - Parallel function application

### Constants

- `STANDARD_BRACKET_TYPES` - Standard bracket types for pseudoknots: `[('(', ')'), ('[', ']'), ('{', '}'), ('<', '>')]`

## Documentation

- **Jupyter Notebooks**: See `notebooks/` directory for detailed examples
  - `01_basic_usage.ipynb` - Basic operations
  - `02_connectivity.ipynb` - Connectivity analysis
  - `03_structure_manipulation.ipynb` - Structure manipulation
  - `04_search_and_analysis.ipynb` - Advanced search
  - `05_json_serialization.ipynb` - JSON operations
  - `06_pandas_integration.ipynb` - Pandas integration
  - `07_parallel_processing.ipynb` - Parallel processing
  - All notebooks have been tested and work with the current version
  - Run `jupyter notebook` from the project root to explore examples
- **API Documentation**: Check docstrings in source code
- **Examples**: All examples in this README are runnable
- **Type Hints**: Full type annotations throughout for better IDE support and type checking

## Development

### Running Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=rna_secstruct --cov-report=html

# Run specific test file
pytest test/test_parser.py

# Run excluding integration tests
pytest -m "not integration"
```

### Code Quality

```bash
# Format code
black rna_secstruct/ test/

# Lint and auto-fix
ruff check rna_secstruct/ test/
ruff check --fix rna_secstruct/ test/

# Type checking
mypy rna_secstruct/

# Run all checks
make check-all
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
  version = {0.1.1},
  year = {2024}
}
```

## Links

- **GitHub**: https://github.com/jyesselm/rna_secstruct
- **Issues**: https://github.com/jyesselm/rna_secstruct/issues
- **Author**: Joe Yesselman (jyesselm@unl.edu)

---

**Note**: This package is designed for non-commercial use. For commercial applications, please contact the author for licensing options.
