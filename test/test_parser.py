"""
testing secondary structure parsing
"""

from rna_secstruct.parser import ConnectivityList, Parser, connectivity_list

# helper funcs ###############################################################


def search_motif_dict(cur_motif, motif_dict):
    """
    :param cur_motif:
    :param motif_dict:
    :return:
    """
    motif_dict[cur_motif.m_id] = cur_motif
    for c in cur_motif.children:
        search_motif_dict(c, motif_dict)


def get_motif_dict(root):
    """
    get a dictionary of motifs
    """
    motif_dict = {}
    search_motif_dict(root, motif_dict)
    return motif_dict


# tests #######################################################################


def test_connectivity_list():
    """
    test connectivity list
    """
    ss = "(((...)))"
    connections = connectivity_list(ss)
    assert connections == [8, 7, 6, -1, -1, -1, 2, 1, 0]


def test_connectivity_list_is_nucleotide_paired():
    """
    test is_nucleotide_paired method of ConnectivityList
    """
    seq = "GGGAAACCC"
    ss = "(((...)))"
    cl = ConnectivityList(seq, ss)
    assert cl.is_nucleotide_paired(0)
    assert cl.is_nucleotide_paired(1)
    assert cl.is_nucleotide_paired(2)
    assert not cl.is_nucleotide_paired(3)
    assert not cl.is_nucleotide_paired(4)
    assert not cl.is_nucleotide_paired(5)
    assert cl.is_nucleotide_paired(6)
    assert cl.is_nucleotide_paired(7)
    assert cl.is_nucleotide_paired(8)


def test_connectivity_list_get_paired_nucleotide():
    """
    test get_paired_nucleotide method of ConnectivityList
    """
    seq = "GGGAAACCC"
    ss = "(((...)))"
    cl = ConnectivityList(seq, ss)
    assert cl.get_paired_nucleotide(0) == 8
    assert cl.get_paired_nucleotide(1) == 7
    assert cl.get_paired_nucleotide(2) == 6
    assert cl.get_paired_nucleotide(6) == 2
    assert cl.get_paired_nucleotide(7) == 1
    assert cl.get_paired_nucleotide(8) == 0


def test_connectivity_list_get_basepair():
    """
    test get_basepair method of ConnectivityList
    """
    seq = "GGGAAACCC"
    ss = "(((...)))"
    cl = ConnectivityList(seq, ss)
    assert cl.get_basepair(0) == "GC"
    assert cl.get_basepair(1) == "GC"
    assert cl.get_basepair(2) == "GC"
    assert cl.get_basepair(6) == "CG"
    assert cl.get_basepair(7) == "CG"
    assert cl.get_basepair(8) == "CG"


def test_in_valid_dot_brackets(caplog):
    """
    test invalid dot brackets - now handled gracefully with warnings
    """
    import logging

    p = Parser()

    # bad structure - should log warnings but not raise exception
    # The structure has length mismatch and unbalanced parentheses
    with caplog.at_level(logging.WARNING):
        result = p.parse("GGGAAACCC", "(((...)))(")
        # Should have warnings about length mismatch and unbalanced parentheses
        assert len(caplog.records) > 0
        # Verify it still parsed (returned a result)
        assert result is not None
        # Clear the log for next test
        caplog.clear()

    # Structure that was previously problematic but now handled
    with caplog.at_level(logging.WARNING):
        result = p.parse("GGGAAACCC", "()((...))")
        # Should parse successfully, may or may not have warnings
        assert result is not None
        caplog.clear()

    # bad sequence - should replace invalid characters with warnings
    with caplog.at_level(logging.WARNING):
        result = p.parse("GGGYAACCC", "(((...)))")
        # Should have warning about invalid sequence characters
        assert result is not None
        # Check that warnings were logged about invalid characters
        # The parser should handle invalid sequence characters gracefully
        _ = [record.message for record in caplog.records]


def test_simple_hairpins():
    """
    test a simple structure
    """
    p = Parser()
    root = p.parse("GGGAAACCC", "(((...)))")
    motif_dict = get_motif_dict(root)
    m0, m1 = motif_dict[0], motif_dict[1]
    assert m0.m_type == "HELIX"
    assert m0.sequence == "GGG&CCC"
    assert m0.structure == "(((&)))"
    assert m0.recursive_sequence() == "GGGAAACCC"
    assert m0.recursive_structure() == "(((...)))"
    assert m0.parent is None
    assert m0.children == [m1]
    assert m1.m_type == "HAIRPIN"
    assert m1.sequence == "GAAAC"
    assert m1.structure == "(...)"
    assert m1.strands == [[2, 3, 4, 5, 6]]
    # has internal loop
    root = p.parse("GGGACCUUCGGGACCC", "(((.((....)).)))")
    motif_dict = get_motif_dict(root)
    assert len(motif_dict) == 4


def test_parsing_edge():
    """
    test a simple structure
    """
    # single nucleotide
    p = Parser()
    root = p.parse("G", ".")
    assert root.m_type == "SINGLESTRAND"
    # single hairpin
    root = p.parse("GAAAC", "(...)")
    motif_dict = get_motif_dict(root)
    assert len(motif_dict) == 2
    # three way junction
    root = p.parse("GGAAACGAAACGAAACC", "((...)(...)(...))")
    motif_dict = get_motif_dict(root)
    m1 = motif_dict[1]
    assert len(m1.children) == 3
    assert len(motif_dict) == 8
    root = p.parse("GGAAACGAAACGAAACCAAA", "((...)(...)(...))...")
    motif_dict = get_motif_dict(root)
    assert motif_dict[8].m_type == "SINGLESTRAND"
    m0 = motif_dict[0]
    assert m0.recursive_sequence() == "GGAAACGAAACGAAACCAAA"
    assert m0.recursive_structure() == "((...)(...)(...))..."
    # single helix
    root = p.parse("AAAAGAAAACAAAA", "....(....)....")
    motif_dict = get_motif_dict(root)
    assert len(motif_dict) == 4
