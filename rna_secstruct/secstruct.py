"""
representation of secondary structure with motif
"""

from rna_secstruct.parser import Parser, is_valid_dot_bracket_str


class SecStruct:
    """
    A class to represent a secondary structure.
    """

    def __init__(self, sequence: str, structure: str):
        """
        :param sequence: a sequence of nucleotides
        :param structure: a dot bracket structure
        """
        self.root = Parser().parse(sequence, structure)
        self.sequence = sequence
        self.structure = structure
        self.motifs = self.__get_motifs(self.root)

    def __add__(self, other):
        """
        add two secondary structures together
        """
        return SecStruct(
            self.sequence + other.sequence, self.structure + other.structure
        )

    def __iter__(self):
        """
        iterate over the motifs
        """
        return iter(self.motifs.values())

    def __getitem__(self, item):
        """
        get a motif by id
        """
        if item not in self.motifs:
            raise ValueError(f"no motif with id {item}")
        return self.motifs[item]

    def __repr__(self):
        return f"{self.sequence}, {self.structure}"

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

    def change_motif(self, m_id, sequence, structure):
        """
        change a motif
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
        full_seq = self.root.recursive_sequence()
        full_ss = self.root.recursive_structure()
        self.sequence = full_seq
        self.structure = full_ss
        self.root = Parser().parse(full_seq, full_ss)
        self.motifs = self.__get_motifs(self.root)

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

    # getters and setters ######################################################
    def get_copy(self):
        """
        get a copy of this secondary structure
        """
        return SecStruct(self.sequence, self.structure)

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
        return len(self.motifs)

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

    def get_motif_by_strand_lengths(self, strand_lengths, m_type=None):
        """
        get a motifs by the length of each strand, strand lengths is an list of n size. Where
        n is the size of each strand. The list is in order of the strand size. For example,
        if the motif is a hairpin it can only have 1 strand, so the list will be [n]. If the
        motif is a junction it can have 2 or more strands, so the list will be [n, m, ...].
        """
        for m in self:
            if m_type is not None and m.type != m_type:
                continue
            # get the length of each strand
            lengths = [len(s) for s in m.strands]
            # compare the lengths
            if lengths == strand_lengths:
                return m
        return None

    def get_twoway_junction_by_topology(self, x_pos, y_pos):
        """
        get a two way junction by topology
        """
        # add 2 two each number in topology
        # to account for the 2 flanking base pairs
        topology = [t + 2 for t in [x_pos, y_pos]]
        return self.get_motif_by_strand_lengths(topology, "JUNCTION")

    def get_motif(self, sequence, structure, min_pos=-1, max_pos=99999):
        """
        get a motif by sequence and structure
        """
        for m in self:
            if m.sequence == sequence and m.structure == structure:
                if m.start_pos >= min_pos and m.end_pos <= max_pos:
                    return m
        return None

    def get_motifs(self, sequence, structure, min_pos=-1, max_pos=99999):
        """
        get motifs by sequence and structure
        """
        motifs = []
        for m in self:
            if m.sequence == sequence and m.structure == structure:
                if m.start_pos >= min_pos and m.end_pos <= max_pos:
                    motifs.append(m)
        return motifs
