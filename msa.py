import sys
import numpy as np
from blosum62 import blosum62
from nuc44 import nuc44


def load_sequences(file_path):
    """Return a list of sequences extracted from the specified file"""
    with open(file_path, "r") as file:
        sequences = [row for row in file.read().splitlines()]
    return sequences


def physico_chemical_group(aminoacid):
    """Return an amino acids's physico-chemical characteristic"""
    hydrophobic = ["A", "V", "I", "L", "M", "F", "Y", "W"]
    electrical = ["R", "H", "K", "D", "E"]
    polar = ["S", "T", "N", "Q"]
    special = ["C", "G", "P"]
    if aminoacid in hydrophobic:
        return "H"
    elif aminoacid in electrical:
        return "E"
    elif aminoacid in polar:
        return "P"
    elif aminoacid in special:
        return "S"
    return "None"


def structure_group(nucleobase):
    """Return the group of a nucleobase based on its structure"""
    purine = ["a", "g"]
    pyrimidine = ["t", "c", "u"]
    if nucleobase in purine:
        return "Pu"
    elif nucleobase in pyrimidine:
        return "Py"
    return "None"


def random_encoding(template):
    """Function to generate a random gap-oriented encoding"""
    encoding = []
    max_seq_length = max([len(seq) for seq in template])

    # Fill the offset in each sequence using random gaps
    for seq in template:
        gene = []
        offset = max_seq_length - len(seq)
        if offset:
            while not len(gene) == offset:
                random_index = np.random.choice(max_seq_length)
                if random_index not in gene:
                    gene.append(random_index)
        encoding.append(sorted(gene))
    return encoding


def decode(encoding, template):
    """Function to decode the encoded information into the original template"""
    alignments = []
    for i in range(len(encoding)):
        alignment = [char for char in template[i]]
        for gap in encoding[i]:
            alignment.insert(gap, "-")
        alignments.append(alignment)
    return np.array(alignments)


def display_msa(alignments, seq_type="protein"):
    """Function to display an alignments"""
    # Make sure that standard python IDLE is used
    try:
        renderer = sys.stdout.shell
    except AttributeError:
        raise RuntimeError("Use IDLE")

    # Compute color code for each column
    column_colorcodes = []
    columns = [alignments[:, i] for i in range(alignments.shape[1])]
    for column in columns:
        if len(set(column)) == 1:
            column_colorcodes.append("hit")
        elif seq_type == "protein":
            if len(set([physico_chemical_group(aa) for aa in column])) == 1:
                column_colorcodes.append("ERROR")
            else:
                column_colorcodes.append("console")
        elif seq_type == "dna":
            if len(set([structure_group(nb) for nb in column])) == 1:
                column_colorcodes.append("ERROR")
            else:
                column_colorcodes.append("console")

    # Display the alignment with respect to the colorcodes
    for seq in alignments:
        for char in range(len(seq)):
            renderer.write(seq[char], column_colorcodes[char])
        print()


def sum_of_pairs_score(encoding, template, seq_type):
    """Return the sum of pairs score of an alignment"""
    alignments = decode(encoding, template)
    # The score of each column is calculated and then summed
    sum_scores = 0
    columns = [alignments[:, i] for i in range(alignments.shape[1])]
    sub_matrices = {"protein": blosum62(), "dna": nuc44()}

    # Match and mismatch scores are calculated using blosum62 and nuc44 matrices
    # and gap score is set to a negative number
    for column in columns:
        score = 0
        for ith_residue in range(len(column)):
            for jth_residue in range(ith_residue + 1, len(column)):
                if not (column[ith_residue] == "-" or column[jth_residue] == "-"):
                    score += sub_matrices[seq_type][(column[ith_residue], column[jth_residue])]
                elif column[ith_residue] == "-" and column[jth_residue] == "-":
                    score += 0
                else:
                    score += -3
        sum_scores += score
    return sum_scores
