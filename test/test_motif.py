"""
tesing the motif module
"""

from rna_secstruct.motif import Motif

def test_simple():
    """
    test a simple structure
    """
    m = Motif("HAIRPIN", [[0, 1, 2, 3, 4]], "GAAAC", "(...)", 0)
    assert m.m_type == "HAIRPIN"
    assert m.sequence == "GAAAC"
    assert m.structure == "(...)"
    assert m.strands == [[0, 1, 2, 3, 4]]
    assert m.recursive_sequence() == "GAAAC"
    assert m.recursive_structure() == "(...)"