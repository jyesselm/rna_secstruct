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

ID: 0, Helix3 GGG&CCC (((&)))<br>
&nbsp;&nbsp;&nbsp;&nbsp;ID: 1, Hairpin3 GAAAC (...)

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

{'m_id': 0, <br>
'm_type': 'HELIX', <br>
'sequence': 'GGG&CCC', <br>
'structure': '(((&)))', <br>
'strands': [[0, 1, 2], [6, 7, 8]], <br>
'positions': [0, 1, 2, 6, 7, 8], <br>
'start_pos': 0, <br>
'end_pos': 8, <br>
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

```python
# recursive functions call all children to get the full sequence and structure of this node and all
# children
m0.recursive_sequence(), m0.recursive_structure()
```

('GGGAAACCC', '(((...)))')

### Working with Secstruct object

```python
struct = SecStruct("GGAAACGAAACGAAACC", "((...)(...)(...))")
```

```python
# getting all helices.
# there are similare get_single_strand(), get_junctions(), get_hairpins()
struct.get_helices()
```

[HELIX,G&C,(&), HELIX,G&C,(&), HELIX,G&C,(&), HELIX,G&C,(&)]

```python
# getting a copy. Always a good idea if you are going to change the structure
# see later sections on design
struct_copy = struct.get_copy()
```

```python
# getting a sub structure
struct = SecStruct("GGGACCUUCGGGACCC", "(((.((....)).)))")
# if I want remove the bottom helix I can simply do
sub_struct = struct.get_sub_structure(1)
# this gets all nodes from 1 and all its children
sub_struct
```

GACCUUCGGGAC, (.((....)).)

#### searching for motifs

```python
# can search for motifs using get_motifs function with search parameters
struct = SecStruct("GGGACCUUCGGGACCC", "(((.((....)).)))")
msp = MotifSearchParams(sequence="GAC&GAC")
struct.get_motifs(msp)
```

[JUNCTION,GAC&GAC,(.(&).)]

```python
msp
```

MotifSearchParams(sequence='GAC&GAC', structure=None, m_type=None, min_pos=0, max_pos=999, min_id=0, max_id=999)

You can search by sequence, structure, type of motifs, the miminal or maximum nucleotide position or the min and max id. These constraints can be useful if you wish to exclude common sequences on the 5' or 3' ends of a sequence of interest

```python
seq = (
    "GGAAGAUCGAGUAGAUCAAAGAGCCUAUGGCUGCCACCCGAGCCCUUGAACUACAGGGAACACUGGAAA"
    "CAGUACCCCCUGCAAGGGCGUUUGACGGUGGCAGCCUAAGGGCUCAAAGAAACAACAACAACAAC"
)
ss = (
    "....((((.....))))...((((((..((((((((((((((((((((.....(((((...((((....)"
    ")))...))))))))))))..)))..))))))))))...))))))...................."
)
struct = SecStruct(seq, ss)
```

```python
msp = MotifSearchParams(m_type="JUNCTION", min_pos=50)
struct.get_motifs(msp)
```

[JUNCTION,GAACA&UACCC,(...(&)...)]

```python
msp = MotifSearchParams(structure="(....)")
struct.get_motifs(msp)

```

[HAIRPIN,GGAAAC,(....)]

```python
# it is also possible to search using motif "tokens" which are the indentifiers generated
# for each motif. For example "Helix4" is any helix of length 4
struct.get_motifs_by_token("Helix4")

```

[HELIX,GAUC&GAUC,((((&)))), HELIX,ACUG&CAGU,((((&))))]

```python
# get a two way junction (Junction2) which has 5 unpaired nucleotides on the first
# strand and 0 on the other
struct.get_motifs_by_token("Junction2_5|0")
```

[JUNCTION,GAACUAC&GC,(.....(&))]

#### changing motifs in secstruct

```python
# using change_motif we can change the sequence or structure of a given motif
# here is a trival example of changing the sequence of the helix which has the id=0
struct = SecStruct("GGGAAACCC", "(((...)))")
struct.change_motif(0, "AGG&CCU", "(((&)))")
struct.sequence, struct.structure
```

('AGGAAACCU', '(((...)))')

```python
# here we can change the loop sequence from a tetraloop to a hexaloop
struct.change_motif(1, "CUUUUUUG", "(......)")
struct.sequence, struct.structure
```

('AGCUUUUUUGCU', '(((......)))')

```python
# can also make larger changes, by replacing a motif with a segment composed
# of many motifs, the object will reparse and create new motifs
struct = SecStruct("GGGAAACCC", "(((...)))")
print("original motifs:")
print(struct.to_str())
print("\nnew motifs")
struct.change_motif(1, "GGGACCUUCGGGACCC", "(((.((....)).)))")
print(struct.to_str())
```

original motifs:
ID: 0, Helix3 GGG&CCC (((&)))
   ID: 1, Hairpin3 GAAAC (...)

new motifs
ID: 0, Helix5 GGGGG&CCCCC (((((&)))))
   ID: 1, Junction2_1|1 GAC&GAC (.(&).)
      ID: 2, Helix2 CC&GG ((&))
          ID: 3, Hairpin4 CUUCGG (....)
