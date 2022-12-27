"""
representation of secondary structure with motif
"""

from rna_secstruct.parser import Parser


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

    def get_sub_structure(self, root_id):
        m = self[root_id]
        seq = m.recursive_sequence()
        struct = m.recursive_structure()
        return SecStruct(seq, struct)

    # getters and setters ######################################################
    def get_num_motifs(self):
        """
        get the number of motifs
        """
        return len(self.motifs)