"""
testing secondary structure parsing
"""
import pytest
from rna_secstruct.parser import Parser, connectivity_list


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


def test_in_valid_dot_brackets():
    """
    test invalid dot brackets
    """
    p = Parser()
    # bad structure
    with pytest.raises(ValueError):
        p.parse("GGGAAACCC", "(((...)))(")
    with pytest.raises(ValueError):
        p.parse("GGGAAACCC", "()((...))")
    # bad sequence
    with pytest.raises(ValueError):
        p.parse("GGGYAACCC", "(((...)))")


def test_simple_hairpin():
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

def test_parsing_edge():
    """
    test a simple structure
    """
    p = Parser()
    root = p.parse("GGG", "...")
    assert root.m_type == "SINGLESTRAND"
