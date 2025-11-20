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


def test_itermotifs():
    struct = SecStruct("GGGAAACCC", "(((...)))")
    inds = []
    for i, motif in struct.itermotifs():
        inds.append(i)
    assert inds == [0, 1]


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


def test_lazy_loading():
    """
    Test that lazy loading works correctly - parsing only happens when motifs are accessed.
    """
    # Create a structure - should not parse immediately
    struct = SecStruct("GGGAAACCC", "(((...)))")
    
    # Accessing sequence/structure should not trigger parsing
    assert struct.sequence == "GGGAAACCC"
    assert struct.structure == "(((...)))"
    
    # __repr__ should not trigger parsing
    repr_str = repr(struct)
    assert "GGGAAACCC" in repr_str
    
    # __add__ should not trigger parsing
    struct2 = SecStruct("AAA", "...")
    combined = struct + struct2
    assert combined.sequence == "GGGAAACCCAAA"
    assert combined.structure == "(((...)))..."
    
    # Accessing motifs should trigger parsing
    # After accessing motifs, they should be cached
    motifs = struct.motifs
    assert len(motifs) == 2
    assert 0 in motifs
    assert 1 in motifs
    
    # Accessing again should use cached version
    motifs2 = struct.motifs
    assert motifs is motifs2  # Should be the same object (cached)
    
    # Accessing _root should also work
    root = struct._root
    assert root is not None
    
    # get_num_motifs should trigger parsing if not already done
    struct3 = SecStruct("GGGAAACCC", "(((...)))")
    num = struct3.get_num_motifs()
    assert num == 2


def test_len():
    """
    Test __len__ method
    """
    struct = SecStruct("GGGAAACCC", "(((...)))")
    assert len(struct) == 9
    assert len(struct) == len(struct.sequence)


def test_slicing():
    """
    Test slicing support in __getitem__
    """
    struct = SecStruct("GGGAAACCC", "(((...)))")
    
    # Test slicing returns new SecStruct
    sliced = struct[2:7]
    assert isinstance(sliced, SecStruct)
    assert sliced.sequence == "GAAAC"
    assert sliced.structure == "(...)"
    
    # Original unchanged
    assert struct.sequence == "GGGAAACCC"
    
    # Test full slice
    full = struct[:]
    assert full.sequence == struct.sequence
    assert full.structure == struct.structure
    assert full is not struct  # Different object
    
    # Test motif access still works
    motif = struct[0]
    assert motif.m_id == 0


def test_split_strands():
    """
    Test split_strands method
    """
    struct = SecStruct("GGG&AAA&CCC", "(((&)))&...")
    strands = struct.split_strands()
    
    assert len(strands) == 3
    assert strands[0].sequence == "GGG"
    assert strands[0].structure == "((("
    assert strands[1].sequence == "AAA"
    assert strands[1].structure == ")))"
    assert strands[2].sequence == "CCC"
    assert strands[2].structure == "..."
    
    # Original unchanged
    assert struct.sequence == "GGG&AAA&CCC"
    
    # Single strand
    struct2 = SecStruct("GGGAAACCC", "(((...)))")
    strands2 = struct2.split_strands()
    assert len(strands2) == 1
    assert strands2[0].sequence == "GGGAAACCC"


def test_insert():
    """
    Test insert method (immutable)
    """
    struct = SecStruct("GGGAAACCC", "(((...)))")
    other = SecStruct("XXX", "...")
    
    # Insert at beginning
    result = struct.insert(0, other)
    assert result.sequence == "XXXGGGAAACCC"
    assert result.structure == "...(((...)))"
    assert struct.sequence == "GGGAAACCC"  # Original unchanged
    
    # Insert in middle
    result2 = struct.insert(3, other)
    assert result2.sequence == "GGGXXXAAACCC"
    assert result2.structure == "(((...)))"  # Structure: "(((" + "..." + "...)))" = "(((...)))" (12 chars)
    
    # Insert at end
    result3 = struct.insert(9, other)
    assert result3.sequence == "GGGAAACCCXXX"
    assert result3.structure == "(((...)))..."
    
    # Test invalid position
    with pytest.raises(ValueError):
        struct.insert(-1, other)
    with pytest.raises(ValueError):
        struct.insert(10, other)


def test_join():
    """
    Test join method (immutable)
    """
    struct1 = SecStruct("GGG", "(((")
    struct2 = SecStruct("AAA", "...")
    
    result = struct1.join(struct2)
    assert result.sequence == "GGG&AAA"
    assert result.structure == "(((&..."
    assert struct1.sequence == "GGG"  # Original unchanged
    assert struct2.sequence == "AAA"  # Original unchanged


def test_replace():
    """
    Test replace method (immutable)
    """
    struct = SecStruct("GGGAAACCC", "(((...)))")
    other = SecStruct("XXX", "...")
    
    # Replace in middle
    result = struct.replace(other, 3)
    assert result.sequence == "GGGXXXCCC"
    assert result.structure == "(((...)))"  # Structure adjusted
    assert struct.sequence == "GGGAAACCC"  # Original unchanged
    
    # Replace at beginning
    result2 = struct.replace(other, 0)
    assert result2.sequence == "XXXAAACCC"
    
    # Test invalid position
    with pytest.raises(ValueError):
        struct.replace(other, -1)
    with pytest.raises(ValueError):
        struct.replace(other, 10)
    with pytest.raises(ValueError):
        struct.replace(SecStruct("XXXXXXXXXX", ".........."), 0)


def test_remove():
    """
    Test remove method (immutable)
    """
    struct = SecStruct("GGGAAACCC", "(((...)))")
    
    # Remove middle region
    result = struct.remove(3, 6)
    assert result.sequence == "GGGCCC"
    assert result.structure == "((()))"  # Adjusted
    assert struct.sequence == "GGGAAACCC"  # Original unchanged
    
    # Remove from beginning
    result2 = struct.remove(0, 3)
    assert result2.sequence == "AAACCC"
    
    # Remove from end
    result3 = struct.remove(6, 9)
    assert result3.sequence == "GGGAAA"
    
    # Test invalid ranges
    with pytest.raises(ValueError):
        struct.remove(-1, 5)
    with pytest.raises(ValueError):
        struct.remove(3, 10)
    with pytest.raises(ValueError):
        struct.remove(5, 3)  # start >= end


def test_subtract():
    """
    Test subtract method (immutable)
    """
    struct = SecStruct("GGGAAACCC", "(((...)))")
    other = SecStruct("AAA", "...")
    
    # Subtract found substructure
    result = struct.subtract(other)
    assert result.sequence == "GGGCCC"
    assert struct.sequence == "GGGAAACCC"  # Original unchanged
    
    # Test not found
    with pytest.raises(ValueError, match="Substructure not found"):
        struct.subtract(SecStruct("XXX", "..."))


def test_to_dict():
    """
    Test to_dict method
    """
    struct = SecStruct("GGGAAACCC", "(((...)))")
    d = struct.to_dict()
    
    assert isinstance(d, dict)
    assert d["sequence"] == "GGGAAACCC"
    assert d["structure"] == "(((...)))"


def test_to_comma_delimited():
    """
    Test to_comma_delimited method
    """
    struct = SecStruct("GGGAAACCC", "(((...)))")
    csv = struct.to_comma_delimited()
    
    assert csv == "GGGAAACCC,(((...)))"
    assert isinstance(csv, str)


def test_immutable_pattern():
    """
    Test that all manipulation methods return new instances (immutable pattern)
    """
    struct = SecStruct("GGGAAACCC", "(((...)))")
    original_id = id(struct)
    
    # Test that operations return new instances
    assert id(struct.insert(3, SecStruct("X", "."))) != original_id
    assert id(struct.join(SecStruct("X", "."))) != original_id
    assert id(struct.replace(SecStruct("X", "."), 0)) != original_id
    assert id(struct.remove(3, 6)) != original_id
    assert id(struct.subtract(SecStruct("AAA", "..."))) != original_id
    assert id(struct[2:7]) != original_id  # Slicing
    
    # Original should be unchanged
    assert struct.sequence == "GGGAAACCC"
    assert id(struct) == original_id
