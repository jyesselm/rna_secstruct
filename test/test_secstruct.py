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
    # trival substructure should be the same as the original
    sub = struct.get_sub_structure(0)
    assert sub.sequence == "GGGACCUUCGGGACCC"
    assert sub[0].sequence == "GGG&CCC"
    assert sub.get_num_motifs() == struct.get_num_motifs()
    sub = struct.get_sub_structure(1)
    assert sub.sequence == "GACCUUCGGGAC"
    sub = struct.get_sub_structure(2)
    assert sub.sequence == "CCUUCGGG"


def test_change_motif_simple():
    """
    testing changing a motif
    """
    struct = SecStruct("GGGAAACCC", "(((...)))")
    struct.change_motif(0, "AGG&CCU", "(((&)))")
    assert struct.sequence == "AGGAAACCU"
    assert struct[0].sequence == "AGG&CCU"
    struct.change_motif(1, "CUUUUUUG", "(......)")
    assert struct.sequence == "AGCUUUUUUGCU"
    assert struct.structure == "(((......)))"
    struct.change_motif(0, "GGG&CCC", "(((&)))")
    assert struct.sequence == "GGGUUUUUUCCC"

def test_change_motif_2():
    """
    testing changing a motif
    """
    struct = SecStruct("GGGACCUUCGGGACCC", "(((.((....)).)))")
    struct.change_motif(1, "GAAAC&GAAAC", "(...(&)...)")
    assert struct.sequence == "GGGAAACCUUCGGGAAACCC"


def test_change_motif_edge_cases():
    """
    testing changing a motif
    """
    struct = SecStruct("GGGACCUUCGGGACCC", "(((.((....)).)))")
    struct.change_motif(1, "GAAAC", "(...)")
    assert struct.sequence == "GGGAAACCC"

    struct = SecStruct("GGAAACGAAACGAAACC", "((...)(...)(...))")
    struct.change_motif(3, "GGGGGAAACCCCC", "(((((...)))))")
    assert struct.sequence == "GGGGGGAAACCCCCGAAACGAAACC"
    assert struct.structure == "((((((...)))))(...)(...))"
    struct = SecStruct("GGAAACGAAACGAAACC", "((...)(...)(...))")
    struct.change_motif(2, "GAGGG&CCACC", "(.(((&)).))")
    assert struct.sequence == "GGAGGGAAACCACCGAAACGAAACC"
    assert struct.structure == "((.(((...)).))(...)(...))"

def test_display():
    """
    testing display
    """
    struct = SecStruct("GGGACCUUCGGGACCC", "(((.((....)).)))")
    print()
    print(struct.to_str())


