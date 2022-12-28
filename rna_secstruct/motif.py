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
        self.__children = []
        self.__end_pos = -1
        self.__parent = None
        self.__positions = []
        self.__m_type = m_type
        self.__m_id = m_id
        self.__start_pos = 9999
        self.__strands = strands
        self.__sequence = sequence
        self.__structure = structure
        self.__token = self.__get_token()
        for strand in strands:
            self.__start_pos = min(min(strand), self.__start_pos)
            self.__end_pos = max(max(strand), self.__end_pos)
            self.__positions.extend(strand)

    def __repr__(self) -> str:
        """
        String representation of just the motif at hand.

        :return: The :class:`str()` representation of the :class:`Motif()`.
        :rtype: :class:`str()`
        """
        return f"{self.__m_type},{self.__sequence},{self.__structure}"

    # setup ####################################################################

    def __get_token(self):
        if self.__m_type == "SINGLESTRAND":
            return f"SingleStrand{len(self.__sequence)}"
        elif self.__m_type == "HAIRPIN":
            return f"Hairpin{len(self.__sequence) - 2}"
        elif self.__m_type == "HELIX":
            spl = self.__structure.split("&")
            return f"Helix{len(spl[0])}"
        elif self.__m_type == "JUNCTION":
            return f"Junction{len(self.__strands)}_" + "|".join(
                [str(len(strand) - 2) for strand in self.__strands]
            )
        raise ValueError(f"Invalid motif type: {self.__m_type}")

    def add_child(self, child) -> None:
        child.__parent = self
        self.__children.append(child)

    def __recursive_build(self, btype):
        result = ""
        if btype == "SEQUENCE":
            value = self.__sequence
        else:
            value = self.__structure
        if self.__m_type == "SINGLESTRAND":
            result = value
            for c in self.__children:
                result += c.__recursive_build(btype)
        elif self.__m_type == "HAIRPIN":
            if self.__parent is not None:
                result = value[1:-1]
            else:
                result = value
        elif self.__m_type == "HELIX":
            result = value
            if len(self.__children) > 0:
                result = value.split("&")
                result = (
                    result[0] + self.__children[0].__recursive_build(btype) + result[1]
                )
            if len(self.__children) > 1:
                result += self.__children[1].__recursive_build(btype)
        elif self.__m_type == "JUNCTION":
            result = ""
            chunks = [subseq[1:-1] for subseq in value.split("&")]
            for chunk, child in zip(chunks, self.__children):
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
    def children(self):
        """
        Returns the children of the motif.
        """
        return self.__children

    @property
    def end_pos(self):
        """
        Returns the end position of the motif.
        """
        return self.__end_pos

    @property
    def parent(self):
        """
        Returns the parent of the motif.
        """
        return self.__parent

    @property
    def m_type(self):
        """
        The type of motif.
        """
        return self.__m_type

    @property
    def m_id(self):
        """
        The motif ID
        """
        return self.__m_id

    @property
    def start_pos(self):
        """
        Returns the start position of the motif.
        """
        return self.__start_pos

    @property
    def strands(self):
        """
        Returns the strands in the motif.
        """
        return self.__strands

    @property
    def sequence(self):
        """
        Returns the sequence of the motif.
        """
        return self.__sequence

    @property
    def structure(self):
        """
        Returns the structure of the motif.
        """
        return self.__structure

    @property
    def token(self):
        """
        Returns the token for the motif.
        """
        return self.__token

    # getters ##################################################################

    def contains(self, position):
        """
        Returns true if the motif contains the given position.
        """
        return position in self.__positions

    def has_parent(self):
        """
        Returns true if the motif has a parent.
        """
        return self.__parent is not None

    def has_children(self):
        """
        Returns true if the motif has children.
        """
        return len(self.__children) > 0

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

    def is_single_strand(self):
        """
        Returns true if the motif is a single strand.
        """
        return self.__m_type == "SINGLESTRAND"

    def num_strands(self):
        """
        Returns the number of strands in the motif.
        """
        return len(self.__strands)

    def to_str(self, depth=0):
        """
        Returns a string representation of the motif.
        """
        id_str = f"ID: {self.__m_id}, "
        pad = "    " * depth
        if not self.has_children():
            return f"{pad}{id_str}{self.__token} {self.__sequence} {self.__structure}"
        else:
            contents = [""]
            for i, child in enumerate(self.__children):
                contents.append(child.to_str(depth + 1))
            children = "\n".join(contents)
        return f"{pad}{id_str}{self.__token} {self.__sequence} {self.__structure}{children}"

    # setters ##################################################################
    @m_type.setter
    def m_type(self, value):
        # check if value is valid
        if value not in ["SINGLESTRAND", "HAIRPIN", "HELIX", "JUNCTION"]:
            raise ValueError(f"Invalid motif type: {value}")
        self.__m_type = value

    @sequence.setter
    def sequence(self, value):
        self.__sequence = value

    @structure.setter
    def structure(self, value):
        self.__structure = value
