# rna_secstruct

[![PYPI status]( https://badge.fury.io/py/rna_secstruct.png)](http://badge.fury.io/py/rna_secstruct)[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

a minimal package for parsing and editing rna secondary structure

## Install

To install rna_secstruct 

```shell
python -m pip install git+https://github.com/jyesselm/rna_secstruct
```


## Features

```python
# load the module, all we need is SecStruct and MotifSearchParams
from rna_secstruct import SecStruct, MotifSearchParams
```

```python
# initiate SecStruct object with sequence and dot bracket notation
struct = SecStruct("GGGAAACCC", "(((...)))")
```

### How to access motifs
Motifs are stored as a dictionary with each given a specific ID. Motifs

```python
struct.motifs
```

    {0: HELIX,GGG&CCC,(((&))), 1: HAIRPIN,GAAAC,(...)}

```python
# can see a formatted version of all motifs with their corresponding ids
print(struct.to_str())
```

    ID: 0, Helix3 GGG&CCC (((&)))
        ID: 1, Hairpin3 GAAAC (...)

```python
# can iterate over motifs
for m in struct:
    print(m)
```

    HELIX,GGG&CCC,(((&)))
    HAIRPIN,GAAAC,(...)

```python
# access first motif which is id `0`
struct[0]
```

    HELIX,GGG&CCC,(((&)))

### Working with motif objects
An overview of all the properties that are stored in a motif

```python
m0 = struct[0]
{
    # what id is the motif
    "m_id" : m0.m_id,
    # the type of motif SINGLESTRAND, HELIX, HAIRPIN, JUNCTION
    "m_type" : m0.m_type,
    # the sequence of the motif, '&' seperates individal strands
    "sequence" : m0.sequence,
    # the structure in dot bracket notatio
    "structure" : m0.structure,
    # the position of each nucleotide in relation the full sequence/structure
    "strands" : m0.strands,
    # all the positions of in a single list
    "positions" : m0.positions,
    # the first position a motif contains
    "start_pos" : m0.start_pos,
    # the last position a motif contains
    "end_pos" : m0.end_pos,
    # the children nodes of this motif
    "children" : m0.children
}
```

    {'m_id': 0,
     'm_type': 'HELIX',
     'sequence': 'GGG&CCC',
     'structure': '(((&)))',
     'strands': [[0, 1, 2], [6, 7, 8]],
     'positions': [0, 1, 2, 6, 7, 8],
     'start_pos': 0,
     'end_pos': 8,
     'children': [HAIRPIN,GAAAC,(...)]}

#### common motif functions

```python
# is this a motif a child of another motif. Can access with m0.parent.
m0.has_parent()
```

    False

```python
# are there nodes under this one. Helices and junctions will almost always have children
m0.has_children()
```

    True

```python
m0.is_single_strand()
```

    False

```python
m0.is_hairpin()
```

    False

```python
m0.is_junction()
```

    False

```python
m0.is_helix()
```

    True

## TODO
