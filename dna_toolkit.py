import collections
from structures import *


def validateSeq(dna_seq):
    """Checks if string is valid DNA sequence"""
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotides:
            return False
    return tmpseq


def countNucFrequency(seq):
    """Counts nucleotide occurence in a sequence"""
    tmpFreqDict = {"A": 0, "C": 0, "G": 0, "T": 0}
    for nuc in seq:
        tmpFreqDict[nuc] += 1
    return tmpFreqDict
    # More pythonic approach
    # return dict(collections.Counter(seq))


def transcription(seq):
    """DNA to RNA translation. Replaces Thymine with Uracil"""
    return seq.replace("T", "U")


def reverse_compliment(seq):
    """Finds complimentary base pairs and reverses list"""
    return "".join([DNA_ReverseCompliment[nuc] for nuc in seq])[::-1]  # Reverses the list
    # More pythonic approach
    # mapping = str.maketrans("ATCG", "TAGC")
    # return seq.translate(mapping)[::-1]


def gc_content(seq):
    """GC content in DNA/RNA sequence"""
    return round((seq.count("C") + seq.count("G")) / len(seq) * 100)


def gc_content_subsec(seq, k = 20):
    """GC Content in a DNA/RNA sub-sequence length k. k = 20 by default"""
    res = []
    for i in range(0, len(seq) - k + 1, k):
        subseq = seq[i:i + k]
        res.append(gc_content(subseq))
    return res