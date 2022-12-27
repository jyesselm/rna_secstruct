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
        self.m_type = m_type
        self.m_id = m_id
        self.strands = strands
        self.sequence = sequence
        self.structure = structure
        self.parent = None
        self.children = []
        self.token = self.__get_token()
        self.depth = -1

    def __repr__(self) -> str:
        """
        String representation of just the motif at hand.

        :return: The :class:`str()` representation of the :class:`Motif()`.
        :rtype: :class:`str()`
        """
        return f"{self.m_type},{self.sequence},{self.structure}"

    # setup ####################################################################

    def __get_token(self):
        if self.m_type == "SINGLESTRAND":
            return f"SingleStrand{len(self.sequence)}"
        elif self.m_type == "HAIRPIN":
            return f"Hairpin{len(self.sequence)}"
        elif self.m_type == "HELIX":
            return f"Helix{len(self.sequence) - 1 // 2}"
        elif self.m_type == "JUNCTION":
            return f"Junction{len(self.strands)}_" + "|".join(
                [str(len(strand) - 2) for strand in self.strands]
            )
        raise ValueError(f"Invalid motif type: {self.m_type}")

    def add_child(self, child) -> None:
        child.parent = self
        self.children.append(child)

    def __recursive_build(self, btype):
        result = ""
        if btype == "SEQUENCE":
            value = self.sequence
        else:
            value = self.structure
        if self.m_type == "SINGLESTRAND":
            result = value
            for c in self.children:
                result += c.__recursive_build(btype)
        elif self.m_type == "HAIRPIN":
            if self.parent is not None:
                result = value[1:-1]
            else:
                result = value
        elif self.m_type == "HELIX":
            result = value.split("&")
            result = result[0] + self.children[0].__recursive_build(btype) + result[1]
            if len(self.children) > 1:
                result += self.children[1].__recursive_build(btype)
        elif self.m_type == "JUNCTION":
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