"""
Represents a motif in an RNA structure.
"""


class Motif:
    """
    A class to represent a motif in a secondary structure.
    """

    # built ins ################################################################

    def __init__(self, m_type, strands, sequence, structure, m_id):
        """
        Setup a new motif object
        :param m_type: type of motif (SINGLESTRAND, HAIRPIN, HELIX, JUNCTION)
        :param strands:
        :param sequence:
        :param structure:
        :param m_id:
        """
        self.__m_type = m_type
        self.m_id = m_id
        self.strands = strands
        self.sequence = sequence
        self.structure = structure
        self.parent = None
        self.children = []
        self.token = self.__get_token()
        self.depth = -1
        self.start_pos = 9999
        self.end_pos = -1
        self.positions = []
        for strand in strands:
            self.start_pos = min(min(strand), self.start_pos)
            self.end_pos = max(max(strand), self.end_pos)
            self.positions.extend(strand)

    def __repr__(self) -> str:
        """
        String representation of just the motif at hand.

        :return: The :class:`str()` representation of the :class:`Motif()`.
        :rtype: :class:`str()`
        """
        return f"{self.__m_type},{self.sequence},{self.structure}"

    # setup ####################################################################

    def __get_token(self):
        if self.__m_type == "SINGLESTRAND":
            return f"SingleStrand{len(self.sequence)}"
        elif self.__m_type == "HAIRPIN":
            return f"Hairpin{len(self.sequence)}"
        elif self.__m_type == "HELIX":
            return f"Helix{len(self.sequence) - 1 // 2}"
        elif self.__m_type == "JUNCTION":
            return f"Junction{len(self.strands)}_" + "|".join(
                [str(len(strand) - 2) for strand in self.strands]
            )
        raise ValueError(f"Invalid motif type: {self.__m_type}")

    def add_child(self, child) -> None:
        child.parent = self
        self.children.append(child)

    def __recursive_build(self, btype):
        result = ""
        if btype == "SEQUENCE":
            value = self.sequence
        else:
            value = self.structure
        if self.__m_type == "SINGLESTRAND":
            result = value
            for c in self.children:
                result += c.__recursive_build(btype)
        elif self.__m_type == "HAIRPIN":
            if self.parent is not None:
                result = value[1:-1]
            else:
                result = value
        elif self.__m_type == "HELIX":
            result = value
            if len(self.children) > 0:
                result = value.split("&")
                result = (
                    result[0] + self.children[0].__recursive_build(btype) + result[1]
                )
            if len(self.children) > 1:
                result += self.children[1].__recursive_build(btype)
        elif self.__m_type == "JUNCTION":
            result = ""
            chunks = [subseq[1:-1] for subseq in value.split("&")]
            for chunk, child in zip(chunks, self.children):
                result += chunk
                result += child.__recursive_build(btype)
            result += chunks[-1]
        return result

    def recursive_sequence(self):
        return self.__recursive_build("SEQUENCE")

    def recursive_structure(self):
        return self.__recursive_build("STRUCTURE")


    # properties ###############################################################
    @property
    def m_type(self):
        """
        The type of motif.
        """
        return self.__m_type


    # getters ##################################################################

    def contains(self, position):
        """
        Returns true if the motif contains the given position.
        """
        return position in self.positions

    def has_parent(self):
        """
        Returns true if the motif has a parent.
        """
        return self.parent is not None

    def has_children(self):
        """
        Returns true if the motif has children.
        """
        return len(self.children) > 0

    def is_junction(self):
        """
        Returns true if the motif is a junction.
        """
        return self.__m_type == "JUNCTION"

    def is_hairpin(self):
        """
        Returns true if the motif is a hairpin.
        """
        return self.__m_type == "HAIRPIN"

    def is_helix(self):
        """
        Returns true if the motif is a helix.
        """
        return self.__m_type == "HELIX"

    def is_singlestrand(self):
        """
        Returns true if the motif is a single strand.
        """
        return self.__m_type == "SINGLESTRAND"


