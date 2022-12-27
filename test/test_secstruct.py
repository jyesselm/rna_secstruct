"""
testing the secondary structure object
"""

from rna_secstruct.secstruct import SecStruct

def test_simple():
    """
    testing a simple case
    """
    struct = SecStruct("GGGAAACCC", "(((...)))")
    assert struct.sequence == "GGGAAACCC"
    assert struct.structure == "(((...)))"
    assert struct.get_num_motifs() == 2
    assert struct[0].sequence == "GGG&CCC"

def test_substruct():
    """
    testing a substructure
    """
    struct = SecStruct("GGGACCUUCGGGACCC", "(((.((....)).)))")
    sub = struct.get_sub_structure(0)
    assert sub.sequence == "GGGACCUUCGGGACCC"
    assert sub[0].sequence == "GGG&CCC"
    assert sub.get_num_motifs() == struct.get_num_motifs()
