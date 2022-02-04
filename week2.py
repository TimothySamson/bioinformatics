import numpy as np
import pandas as pd
import itertools
from main import ReverseComplement

def HammingDistance(p, q):
    if len(p) != len(q):
        raise ValueError("length of p equals length of q")

    count = 0
    n = len(p)
    for i in range(0, n):
        if p[i] != q[i]:
            count += 1

    return count

# calculates #G - #C from the first i places in the genome
def SkewArray(genome, i=0):
    if i == 0:
        i = len(genome)

    count = 0
    skews = [0]
    for k in range(0, i):
        if genome[k].lower() == "c":
            count -= 1
        elif genome[k].lower() == "g":
            count += 1
        skews.append(count)

    return skews

def MinSkew(genome):
    skewArray = np.array(SkewArray(genome))
    minSkew = min(skewArray)
    return np.where(skewArray == minSkew)[0]

def ApproximatePatternIndices(pattern, genome, d):
    n = len(genome)
    k = len(pattern)
    indices = []
    for i in range(0, n - k + 1):
        subgenome = genome[i: i+k]
        if HammingDistance(pattern, subgenome) <= d:
            indices.append(i)

    return indices

def Count(text, pattern, d):
    n = len(text)
    k = len(pattern)
    count = 0
    for i in range(n-k+1):
        if HammingDistance(text[i: i+k], pattern) <= d:
            count += 1

    return count

# Returns all versions of `text` with hamming distances at most `d`
def HammingNeighbors(text, d):
    n = len(text)
    neighbors = set()

    if d == 0:
        return {text}
    if n == 1:
        return {"A", "T", "G", "C"}

    suffix = text[1:]
    suffixNeighbors = HammingNeighbors(suffix, d)
    for suffixNeighbor in suffixNeighbors:
        if HammingDistance(suffixNeighbor, suffix) < d:
            for nucleotide in "ATGC":
                neighbors.add(nucleotide + suffixNeighbor)
        else:
            neighbors.add(text[0] + suffixNeighbor)

    return neighbors

def HammingAndRCVariants(text, d):
    variants = list(HammingNeighbors(text, d))
    variants.extend([ReverseComplement(x) for x in variants])
    return variants

# The mismatch function takes in a subtext and returns a neighborhood of text
def FrequentWordsWithMismatches(text, k, mismatch):
    freqMap = {}
    n = len(text)
    subtexts = [text[i: i + k] for i in range(n - k + 1)]

    for subtext in subtexts:
        for variant in mismatch(subtext):
            if variant in freqMap:
                freqMap[variant] += 1
            else:
                freqMap[variant] = 1

    freqMap = pd.Series(freqMap)
    return freqMap

if __name__ == "__main__":
    with open("dataset_9_10.txt") as file:
        text = file.readline().strip()
        k, d = map(int, file.readline().strip().split())

        print(FrequentWordsWithMismatches(text, k, lambda x: HammingAndRCVariants(x, d)).sort_values())


    # with open("input_5.txt") as file:
    #     text = file.readline().strip()
    #     k, d = map(int, file.readline().strip().split())
    #
    #     print(FrequentWordsWithMismatches(text, k, d).sort_values())

    # text = "ATGC" # Reverse complement: GCAT






