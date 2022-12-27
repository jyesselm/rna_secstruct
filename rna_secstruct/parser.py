"""
A simple parser of rna secondary structure inspired by `rna_library` code written
by Chris Jurich
"""
import re
from typing import List


def is_valid_dot_bracket_str(structure: str) -> bool:
    """
    Checks if a structure is a valid dot-bracket structure containing only
    '(', '.' or ')' characters.
    :param: str structure: dot bracket structure
    """
    lparen_ct = 0
    for ch in structure:
        if ch == "(":
            lparen_ct += 1
        elif ch == ")":
            lparen_ct -= 1
        elif ch == ".":
            continue
        else:
            raise ValueError(
                f"{ch} is invalid in a dot-bracket structure. Only '(', '.' "
                f"and ')' are allowed"
            )

        if lparen_ct < 0:
            raise ValueError(f"{structure} is an unbalanced structure")

    if lparen_ct != 0:
        raise ValueError(f"{structure} is an unbalanced structure")

    for ii in range(3):
        invalid = "(" + "." * ii + ")"
        if structure.find(invalid) != -1:
            raise ValueError(f"{structure} has a hairpin that is too small")

    return True


def connectivity_list(structure: str) -> List[int]:
    """
    Generates a connectivity list or pairmap from a dot-bracket secondary structure.
    The list has the index of a positions complement, is it is a '.' it will have
    a -1 instead.

    :param structure: a dot-bracket structure
    :rtype: list[int]
    :raises TypeError: if the number of left parentheses exceeds the number of
     right parentheses
    """
    connections, pairs = [-1] * len(structure), []
    for index, db in enumerate(structure):
        if db == "(":
            pairs.append(index)
        elif db == ")":
            complement = pairs.pop()
            connections[complement] = index
            connections[index] = complement
    if len(pairs):
        raise TypeError("Unbalanced parentheses in structure")
    return connections


def is_circular(start, connections):
    it = start + 1
    while True:
        while it < len(connections) and connections[it] == -1:
            it += 1
        if it == len(connections):
            return False

        it = connections[it] + 1
        if it == start or it < start:
            return True


class Motif:
    """
    A class to represent a motif in a secondary structure.
    """

    def __init__(self, m_type, strands, sequence, structure, m_id):
        self.m_type = str(m_type)
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

    def recursive_sequence(self):
        result = ""
        if self.m_type == "SINGLESTRAND":
            result = self.sequence
            for c in self.children:
                result += c.recursive_sequence()
        elif self.m_type == "HAIRPIN":
            result = self.sequence[1:-1]
        elif self.m_type == "HELIX":
            result = self.sequence.split("&")
            result = result[0] + self.children[0].recursive_sequence() + result[1]
            if len(self.children) > 1:
                result += self.children[1].recursive_sequence()
        elif self.m_type == "JUNCTION":
            result = ""
            seq_chunks = [subseq[1:-1] for subseq in self.sequence.split("&")]
            for seq_chunk, child in zip(seq_chunks, self.children):
                result += seq_chunk
                result += child.recursive_sequence()
            result += seq_chunks[-1]
        return result

    def recursive_structure(self):
        result = ""
        if self.m_type == "SINGLESTRAND":
            result = self.structure
            for c in self.children:
                result += c.recursive_structure()
        elif self.m_type == "HAIRPIN":
            result = self.structure[1:-1]
        elif self.m_type == "HELIX":
            result = self.structure.split("&")
            result = result[0] + self.children[0].recursive_structure() + result[1]
            if len(self.children) > 1:
                result += self.children[1].recursive_structure()
        elif self.m_type == "JUNCTION":
            result = ""
            secstruct_chunks = [subseq[1:-1] for subseq in self.structure.split("&")]
            for secstruct_chunk, child in zip(secstruct_chunks, self.children):
                result += secstruct_chunk
                result += child.recursive_structure()
            result += secstruct_chunks[-1]
        return result


class Parser:
    """
    A class to parse secondary structure into motifs
    """

    def __init__(self):
        self.motif_id = 0

    def parse(self, sequence, structure):
        self.motif_id = 0
        self.__check_to_see_if_inputs_valid(sequence, structure)
        connections = connectivity_list(structure)
        return self.__get_motifs(sequence, structure, connections, 0)

    def __check_to_see_if_inputs_valid(self, sequence, structure):
        """
        check to see if the inputs are valid
        :param sequence: a sequence of nucleotides
        :param structure:
        :return: None
        """
        # enforce upper case sequence and is RNA
        sequence = sequence.upper().replace("T", "U")
        if len(sequence) != len(structure):
            raise ValueError(
                f"sequence and structure are not the same length:"
                f" {sequence} {structure}"
            )
        if not re.match(r"^[ACGUTN]+$", sequence):
            raise ValueError(f"sequence contains invalid characters: {sequence}")
        if not re.match(r"^[().]+$", structure):
            raise ValueError(f"structure contains invalid characters: {structure}")
        is_valid_dot_bracket_str(structure)

    def __get_motifs(self, sequence, structure, connections, start):
        """
        get the motifs from the structure
        :param sequence: a sequence of nucleotides
        :param structure: a dot bracket structure
        :param connections: a list of connections
        :return: a list of motifs
        """
        motifs = []
        if start >= len(connections):
            return None
        if connections[start] == -1:
            return self.__get_single_strand(sequence, structure, connections, start)
        return self.__get_helix(sequence, structure, connections, start)

    def __get_single_strand(self, sequence, structure, connections, start):
        """
        get a single strand
        :param sequence: a sequence of nucleotides
        :param connections: a list of connections
        :param start: the start of the single strand
        :return: a single strand
        """
        # how many nucleotides are in the single strand
        single_strand_count = 0
        while (
            start + single_strand_count < len(connections)
            and connections[start + single_strand_count] == -1
        ):
            single_strand_count += 1
        strand = list(range(start, start + single_strand_count))
        sstrand = Motif(
            "SINGLESTRAND",
            strand,
            sequence[start : start + single_strand_count],
            structure[start : start + single_strand_count],
            self.motif_id,
        )
        self.motif_id += 1
        if start + single_strand_count < len(connections):
            sstrand.add_child(
                self.__get_motifs(
                    sequence, structure, connections, start + single_strand_count
                )
            )
        return sstrand

    def __get_helix(self, sequence, structure, connections, start):
        """
        get a helix or junction
        :param sequence: a sequence of nucleotides
        :param connections: a list of connections
        :param start: the start of the helix or junction
        :return: a helix or junction
        """
        helix_len = self.__get_helix_length(connections, start)
        lhs, rhs = [], []
        for index in range(start, start + helix_len):
            lhs.append(index)
            rhs.append(connections[index])
        rhs.reverse()
        seq1, ss1 = self.__get_seq_and_ss_from_strand(sequence, structure, lhs)
        seq2, ss2 = self.__get_seq_and_ss_from_strand(sequence, structure, rhs)
        helix = Motif(
            "HELIX", [lhs, rhs], f"{seq1}&{seq2}", f"{ss1}&{ss2}", self.motif_id
        )
        self.motif_id += 1
        if connections[start + helix_len - 1] > start:
            helix.add_child(
                self.__get_junction_or_hairpin(
                    sequence, structure, connections, start + helix_len - 1
                )
            )

        if not is_circular(rhs[-1], connections):
            motif = self.__get_motifs(sequence, structure, connections, rhs[-1] + 1)
            if motif is not None:
                helix.add_child(motif)
        return helix

    def __get_junction_or_hairpin(self, sequence, structure, connections, start):
        """
        get a junction or hairpin
        :param sequence:
        :param structure:
        :param connections:
        :param start:
        :return:
        """
        strands = []
        pos = start
        # pos should be the first opening pair of a junction or hairpin
        if structure[pos] != "(":
            raise ValueError(f"expected ( at position {pos}")
        while True:
            next_strand = [pos]
            pos += 1
            while connections[pos] == -1:
                next_strand.append(pos)
                pos += 1
            next_strand.append(pos)
            strands.append(next_strand)
            pos = connections[pos]
            # made a complete circle
            if pos == start:
                break
        # is a junction
        self.motif_id += 1
        if len(strands) > 1:
            seq_and_ss = [
                self.__get_seq_and_ss_from_strand(sequence, structure, strand)
                for strand in strands
            ]
            seq = "&".join([seq for seq, ss in seq_and_ss])
            ss = "&".join([ss for seq, ss in seq_and_ss])
            m = Motif("JUNCTION", strands, seq, ss, self.motif_id - 1)
            for strand in strands[:-1]:
                m.add_child(
                    self.__get_motifs(sequence, structure, connections, strand[-1])
                )
            return m
        else:
            seq, ss = self.__get_seq_and_ss_from_strand(sequence, structure, strands[0])
            return Motif("HAIRPIN", strands, seq, ss, self.motif_id - 1)

    def __get_helix_length(self, connections, start):
        """
        get the length of a helix
        :param connections: a list of connections
        :param start: the start of the helix
        :return: the length of the helix
        """
        complement = connections[start]
        length = 0
        while (
            connections[start + length] == complement - length
            and connections[complement - length] == start + length
        ):
            length += 1
        return length

    def __get_seq_and_ss_from_strand(self, sequence, structure, strand):
        """
        get the sequence and secondary structure from a strand
        :return: the sequence and secondary structure
        """
        seq = "".join([sequence[i] for i in strand])
        ss = "".join([structure[i] for i in strand])
        return seq, ss
