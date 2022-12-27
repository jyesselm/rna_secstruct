"""
tesing the motif module
"""

from rna_secstruct.motif import Motif


def test_simple():
    """
    test a simple structure
    """
    m = Motif("SINGLESTRAND", [[0, 1, 2, 3, 4]], "AAAAA", ".....", 0)
    assert m.recursive_sequence() == "AAAAA"
    assert m.recursive_structure() == "....."
    assert m.is_singlestrand()
    assert not m.is_helix()
    assert not m.is_hairpin()
    assert not m.is_junction()
    # hairpin test
    m = Motif("HAIRPIN", [[0, 1, 2, 3, 4]], "GAAAC", "(...)", 0)
    assert m.m_type == "HAIRPIN"
    assert m.sequence == "GAAAC"
    assert m.structure == "(...)"
    assert m.strands == [[0, 1, 2, 3, 4]]
    assert m.recursive_sequence() == "GAAAC"
    assert m.recursive_structure() == "(...)"
    # helix test
    m = Motif("HELIX", [[0, 1, 2], [3, 4, 5]], "GGG&CCC", "(((&)))", 0)
    assert m.recursive_sequence() == "GGG&CCC"
    assert m.recursive_structure() == "(((&)))"
    # junction test
    m = Motif("JUNCTION", [[0, 1, 2], [3, 4, 5]], "GAG&CAC", "(.(&).)", 0)
    # TODO fix this
    # print(m.recursive_sequence())


def test_other_getters():
    """
    test other getting functions
    """
    m = Motif("SINGLESTRAND", [[0, 1, 2, 3, 4]], "AAAAA", ".....", 0)
    assert m.end_pos == 4
    assert m.start_pos == 0
    assert m.contains(0)
    assert not m.has_parent()
    assert not m.has_children()
    m = Motif("JUNCTION", [[0, 1, 2], [3, 4, 5]], "GAG&CAC", "(.(&).)", 0)
    assert m.end_pos == 5
    assert m.start_pos == 0
    assert m.contains(0)


def test_simple_build_up():
    """
    test a simple structure
    """
    m1 = Motif("SINGLESTRAND", [[0, 1, 2]], "AAA", "...", 0)
    m2 = Motif("HELIX", [[0, 1, 2], [3, 4, 5]], "GGG&CCC", "(((&)))", 0)
    m3 = Motif("HAIRPIN", [[0, 1, 2, 3, 4]], "GAAAC", "(...)", 0)
    m1.add_child(m2)
    m2.add_child(m3)
    assert m1.recursive_sequence() == "AAAGGGAAACCC"
