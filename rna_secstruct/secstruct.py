"""
representation of secondary structure with motif
"""

from typing import List, Optional
from dataclasses import dataclass
from rna_secstruct.parser import Parser, is_valid_dot_bracket_str
from rna_secstruct.motif import Motif


@dataclass(order=True)
class MotifSearchParams:
    """
    params for rna design
    """

    sequence: Optional[str] = None
    structure: Optional[str] = None
    m_type: Optional[str] = None
    min_pos: int = 0
    max_pos: int = 999
    min_id: int = 0
    max_id: int = 999


class SecStruct:
    """
    A class to represent a secondary structure.
    """

    def __init__(self, sequence: str, structure: str):
        """
        :param sequence: a sequence of nucleotides
        :param structure: a dot bracket structure
        """
        self.__root = Parser().parse(sequence, structure)
        self.__motifs = self.__get_motifs(self.__root)
        self.__sequence = sequence
        self.__structure = structure

    def __add__(self, other):
        """
        add two secondary structures together
        """
        return SecStruct(
            self.__sequence + other.__sequence, self.__structure + other.__structure
        )

    def itermotifs(self):
        """
        iterate over motifs
        """
        return iter(self.motifs.items())

    def __iter__(self):
        """
        iterate over the motifs
        """
        return iter(self.motifs.values())

    def __getitem__(self, item) -> Motif:
        """
        get a motif by id
        """
        if item not in self.motifs:
            raise ValueError(f"no motif with id {item}")
        return self.motifs[item]

    def __repr__(self):
        return f"{self.__sequence}, {self.__structure}"

    def __get_motifs(self, root):
        """
        get a list of motifs
        """
        motifs = {}
        self.__collect_motifs(root, motifs)
        return motifs

    def __collect_motifs(self, cur_motif, motifs):
        """
        collect motifs
        """
        motifs[cur_motif.m_id] = cur_motif
        for c in cur_motif.children:
            self.__collect_motifs(c, motifs)

    # design ###################################################################

    def change_motif_sequence(self, m_id, seqeuence):
        """
        change the sequence of a motif
        :param m_id:
        :param seqeuence:
        :return:
        """
        pass

    def change_motif(self, m_id, sequence, structure) -> None:
        """
        change a motif sequence and secondary structure triggering a reparse
        of secondary structure.
        :param m_id: the id of the motif to change
        :param sequence: the new sequence
        :param structure: the new structure
        :return: None
        """
        is_valid_dot_bracket_str(structure)
        if len(sequence) != len(structure):
            raise ValueError("sequence and structure must be the same length")
        if sequence.count("&") != structure.count("&"):
            raise ValueError(
                "sequence and structure must have the same number of strands"
            )
        m = self[m_id]
        # logic to change motif type only can in number of strands
        if m.num_strands() < sequence.count("&") + 1:
            raise ValueError("cannot add strands to a motif")
        m_type = m.m_type
        if m.num_strands() > sequence.count("&") + 1 and sequence.count("&") == 0:
            if structure.count("(") == 0:
                m_type = "SINGLESTRAND"
            else:
                m_type = "HAIRPIN"
        if m.has_parent() and m.parent.is_helix():
            self.__change_inner_flanking(m.parent, sequence[0] + sequence[-1])
        strands = sequence.split("&")
        for s1, s2, c in zip(strands[:-1], strands[1:], m.children):
            self.__change_outer_flanking(c, s1[-1] + s2[0])
        m.sequence = sequence
        m.structure = structure
        m.m_type = m_type
        full_seq = self.__root.recursive_sequence()
        full_ss = self.__root.recursive_structure()
        self.__sequence = full_seq
        self.__structure = full_ss
        self.__root = Parser().parse(full_seq, full_ss)
        self.__motifs = self.__get_motifs(self.__root)

    def __change_inner_flanking(self, m, new_cp: str):
        """
        change the inner flanking sequence of a motif
        """
        it = m.sequence.find("&")
        tks = list(m.sequence)
        tks[it - 1], tks[it + 1] = new_cp[0], new_cp[1]
        m.sequence = "".join(tks)

    def __change_outer_flanking(self, m, new_cp: str):
        """
        change the outer flanking sequence of a motif
        """
        tks = list(m.sequence)
        tks[0], tks[-1] = new_cp[0], new_cp[1]
        m.sequence = "".join(tks)

    # search ###################################################################

    def __get_motifs_by_params(self, msp: MotifSearchParams) -> List[Motif]:
        """
        get a list of motifs by params
        :param msp: parameters
        :return:
        """
        motifs = []
        for m in self:
            if msp.m_type is not None and m.m_type != msp.m_type:
                continue
            if msp.sequence is not None and msp.sequence != m.sequence:
                continue
            if msp.structure is not None and msp.structure != m.structure:
                continue
            if m.m_id < msp.min_id or m.m_id > msp.max_id:
                continue
            if m.start_pos < msp.min_pos or m.end_pos > msp.max_pos:
                continue
            motifs.append(m)
        return motifs

    def get_motifs_by_strand_lengths(
        self, strand_lengths, msp: Optional[MotifSearchParams] = None
    ) -> List[Motif]:
        """
        get a motifs by the length of each strand, strand lengths is an list of
        n size. Where n is the size of each strand. The list is in order of the
        strand size. For example, if the motif is a hairpin it can only have
        1 strand, so the list will be [n]. If the motif is a junction it can
        have 2 or more strands, so the list will be [n, m, ...].

        :param strand_lengths: the list of strand lengths, e.g. [3, 4, 5]
        :param msp: optional motif search params
        :return: a list of motifs
        """
        if msp is None:
            msp = MotifSearchParams()
        selected_motifs = self.__get_motifs_by_params(msp)
        motifs = []
        for m in selected_motifs:
            # get the length of each strand
            lengths = [len(s) for s in m.strands]
            # compare the lengths
            if lengths == strand_lengths:
                motifs.append(m)
        return motifs

    def get_twoway_junctions_by_topology(
        self, x_pos, y_pos, msp: Optional[MotifSearchParams] = None
    ) -> List[Motif]:
        """
        get a two way junction by topology
        """
        # add 2 two each number in topology
        # to account for the 2 flanking base pairs
        if msp is None:
            msp = MotifSearchParams(m_type="JUNCTION")
        else:
            if msp.m_type is not None and msp.m_type != "JUNCTION":
                raise ValueError("m_type must be JUNCTION")
            msp.m_type = "JUNCTION"

        topology = [t + 2 for t in [x_pos, y_pos]]
        return self.get_motifs_by_strand_lengths(topology, msp)

    def get_motifs(self, msp) -> List[Motif]:
        """
        get motifs by sequence and structure
        """
        return self.__get_motifs_by_params(msp)

    def get_motifs_by_token(
        self, token, msp: Optional[MotifSearchParams] = None
    ) -> List[Motif]:
        """
        get a list of motifs by a token
        :param token: the token to search for
        :param msp: optional motif search params
        :return: a list of motifs
        """
        if msp is None:
            msp = MotifSearchParams()
        selected_motifs = self.__get_motifs_by_params(msp)
        motifs = []
        for m in selected_motifs:
            if token == m.token:
                motifs.append(m)
        return motifs

    # properites ###############################################################
    @property
    def motifs(self):
        """
        get motifs
        """
        return self.__motifs

    @property
    def sequence(self):
        """
        get the sequence
        """
        return self.__sequence

    @property
    def structure(self):
        """
        get the structure
        """
        return self.__structure

    # getters and setters ######################################################
    def get_copy(self):
        """
        get a copy of this secondary structure
        """
        return SecStruct(self.__sequence, self.__structure)

    def get_hairpins(self):
        """
        get hairpin motifs
        """
        return [m for m in self if m.is_hairpin()]

    def get_helices(self):
        """
        get helix motifs
        """
        return [m for m in self if m.is_helix()]

    def get_junctions(self):
        """
        get junctions
        """
        return [m for m in self if m.is_junction()]

    def get_single_strands(self):
        """
        get single strands
        """
        return [m for m in self if m.is_single_strand()]

    def get_num_motifs(self):
        """
        get the number of motifs
        """
        return len(self.__motifs)

    def get_sub_structure(self, root_id):
        m = self[root_id]
        seq = m.recursive_sequence()
        struct = m.recursive_structure()
        # fix dropping the first basepair in the structure
        if m.has_parent() and m.parent.is_helix():
            seq_spl = m.parent.sequence.split("&")
            seq = seq_spl[0][-1] + seq + seq_spl[1][0]
            struct = "(" + struct + ")"
        return SecStruct(seq, struct)

    def to_str(self):
        """
        get a string representation of this secondary structure
        """
        return self.__root.to_str()
