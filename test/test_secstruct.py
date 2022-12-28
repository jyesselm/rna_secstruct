"""
testing the secondary structure object
"""

import os
import pytest
from rna_secstruct.secstruct import SecStruct, MotifSearchParams

# get the current directory
CUR_DIR = os.path.dirname(os.path.realpath(__file__))


def test_simple():
    """
    testing a simple case
    """
    struct = SecStruct("GGGAAACCC", "(((...)))")
    assert struct.sequence == "GGGAAACCC"
    assert struct.structure == "(((...)))"
    assert struct.get_num_motifs() == 2
    assert struct[0].sequence == "GGG&CCC"


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


def test_copy():
    """
    test deep copy
    """
    struct = SecStruct("GGGACCUUCGGGACCC", "(((.((....)).)))")
    struct_copy = struct.get_copy()
    assert struct_copy.sequence == struct.sequence
    assert struct_copy.structure == struct.structure
    assert struct_copy.get_num_motifs() == struct.get_num_motifs()
    # test that changing the copy does not change the original
    struct_copy.change_motif(0, "AGG&CCU", "(((&)))")
    assert struct_copy.sequence != struct.sequence


def test_get_hairpins():
    """
    testing get hairpins
    """
    struct = SecStruct("GGAAACGAAACGAAACC", "((...)(...)(...))")
    hps = struct.get_hairpins()
    assert len(hps) == 3
    assert hps[0].sequence == "GAAAC"


def test_get_helices():
    """
    testing get helices
    """
    struct = SecStruct("GGAAACGAAACGAAACC", "((...)(...)(...))")
    helices = struct.get_helices()
    assert len(helices) == 4


def test_get_junctions():
    """
    testing getting junctions
    """
    struct = SecStruct("GGAAACGAAACGAAACC", "((...)(...)(...))")
    juncs = struct.get_junctions()
    assert len(juncs) == 1


def test_get_single_strands():
    """
    testing getting single strands
    """
    struct = SecStruct("AAAGGAAACGAAACGAAACCAAA", "...((...)(...)(...))...")
    sstrands = struct.get_single_strands()
    assert len(sstrands) == 2


def test_sub_structure():
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


def test_get_motif_by_strand_lengths():
    """
    testing getting a motif by strand lengths
    """
    struct = SecStruct("GGGACCUUCGGGACCC", "(((.((....)).)))")
    motifs = struct.get_motifs_by_strand_lengths([3, 3])
    assert len(motifs) == 2
    msp = MotifSearchParams(m_type="HELIX")
    motifs = struct.get_motifs_by_strand_lengths([3, 3], msp)
    assert len(motifs) == 1
    msp = MotifSearchParams(min_pos=1)
    motifs = struct.get_motifs_by_strand_lengths([3, 3], msp)
    assert len(motifs) == 1


def test_get_twoway_junctions():
    """
    testing getting two way junctions
    """
    struct = SecStruct("GGGACCUUCGGGACCC", "(((.((....)).)))")
    motifs = struct.get_twoway_junctions_by_topology(1, 1)
    assert len(motifs) == 1
    # cannot supply a m_type other than JUNCTION or nothing at all
    # this seems like a problem or bad design choice
    msp = MotifSearchParams(m_type="HELIX")
    with pytest.raises(ValueError):
        struct.get_twoway_junctions_by_topology(1, 1, msp)


def test_get_motifs():
    """
    test get motifs
    """
    # example from mttr-6-alt-h1 (C000T)
    seq = (
        "GGAAGAUCGAGUAGAUCAAAGAGCCUAUGGCUGCCACCCGAGCCCUUGAACUACAGGGAACACUGGAAA"
        "CAGUACCCCCUGCAAGGGCGUUUGACGGUGGCAGCCUAAGGGCUCAAAGAAACAACAACAACAAC"
    )
    ss = (
        "....((((.....))))...((((((..((((((((((((((((((((.....(((((...((((....)"
        ")))...))))))))))))..)))..))))))))))...))))))...................."
    )
    struct = SecStruct(seq, ss)
    msp = MotifSearchParams(m_type="HELIX")
    motifs = struct.get_motifs(msp)
    assert len(motifs) == 7
    msp = MotifSearchParams(sequence="GAACA&UACCC")
    motifs = struct.get_motifs(msp)
    assert len(motifs) == 1
    msp = MotifSearchParams(m_type="JUNCTION", min_pos=50)
    motifs = struct.get_motifs(msp)
    assert len(motifs) == 1
    msp = MotifSearchParams(structure="(....)")
    motifs = struct.get_motifs(msp)
    assert len(motifs) == 1


def test_get_motifs_by_token():
    """
    testing getting a motifs by token
    """
    # example from mttr-6-alt-h1 (C000T)
    seq = (
        "GGAAGAUCGAGUAGAUCAAAGAGCCUAUGGCUGCCACCCGAGCCCUUGAACUACAGGGAACACUGGAAA"
        "CAGUACCCCCUGCAAGGGCGUUUGACGGUGGCAGCCUAAGGGCUCAAAGAAACAACAACAACAAC"
    )
    ss = (
        "....((((.....))))...((((((..((((((((((((((((((((.....(((((...((((....)"
        ")))...))))))))))))..)))..))))))))))...))))))...................."
    )
    struct = SecStruct(seq, ss)
    motifs = struct.get_motifs_by_token("Helix4")
    assert len(motifs) == 2
    motifs = struct.get_motifs_by_token("Junction2_5|0")
    assert len(motifs) == 1
    motifs = struct.get_motifs_by_token("Junction2_0|2")
    assert len(motifs) == 2


def test_display():
    """
    testing display
    """
    struct = SecStruct("GGGACCUUCGGGACCC", "(((.((....)).)))")
    struct_str = struct.to_str()
    with open(os.path.join(CUR_DIR, "resources", "test_display.txt")) as fin:
        assert struct_str == fin.read().strip()
