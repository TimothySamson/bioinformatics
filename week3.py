import numpy as np
import pandas as pd
import math
from week2 import HammingNeighbors, HammingDistance


def MotifEnumeration(dnas, k, d):
    n = len(dnas[0])
    motifs = set()
    for i in range(n-k+1):
        for nbr in HammingNeighbors(dnas[0][i: i+k], d):
            found = False
            for dna in dnas[1:]:
                found = False
                for j in range(n-k+1):
                    subdna = dna[j:j+k]
                    if HammingDistance(nbr, subdna) <= d:
                        found = True
                        break
                if not found:
                    break
            if found:
                motifs.add(nbr)
    return motifs

def Consensus(dnas):
    if not isinstance(dnas, pd.DataFrame):
        dnas = dna_to_dataframe(dnas)
    return dnas.mode().iloc[0]


def dna_to_dataframe(dnas):
    return pd.DataFrame(map(list, dnas))

# PANDAS TOO SLOW T_T

# def Score(dnas):
#     if not isinstance(dnas, pd.DataFrame):
#         dnas = dna_to_dataframe(dnas)
#     return (Consensus(dnas) != dnas).to_numpy().sum()


# def Count(dnas):
#     if not isinstance(dnas, pd.DataFrame):
#         dnas = dna_to_dataframe(dnas)
#
#     countMatrix = dnas.apply(pd.value_counts).fillna(0)
#     for nucleo in "ATGC":
#         if nucleo not in countMatrix.index:
#             countMatrix.loc[nucleo] = 0
#
#     return countMatrix


# def Profile(dnas, pseudo=False):
#     if not isinstance(dnas, pd.DataFrame):
#         dnas = dna_to_dataframe(dnas)
#
#     count = Count(dnas)
#     if pseudo:
#         count = count.applymap(lambda x: x + 1)
#     return count / count[0].sum()


def Entropy(dnas):
    if not isinstance(dnas, pd.DataFrame):
        dnas = dna_to_dataframe(dnas)

    dnas = Profile(dnas).applymap(lambda x: -x * np.log2(x) if x > 0 else 0)

    return dnas.sum().sum()


def distance(pattern, text):
    n = len(text)
    k = len(pattern)

    if n < k:
        raise ValueError("pattern has to have less length than text")

    return min([HammingDistance(pattern, text[i: i + k]) for i in range(n - k + 1)])


def kmers(k):
    if k == 1:
        return list("ATGC")

    suffixes = kmers(k - 1)
    res = []
    for nucleo in "ATGC":
        res.extend([nucleo + suffix for suffix in suffixes])

    return res


def MedianString(dnas, k):
    bestDist = math.inf
    bestMotifs = []

    for kmer in kmers(k):
        dist = sum([distance(kmer, dna) for dna in dnas])
        if dist < bestDist:
            bestDist = dist
            bestMotifs = [kmer]
        if dist == bestDist:
            bestMotifs.append(kmer)

    return bestMotifs

# def Pr(motif, profile):
#     if len(profile.columns) != len(motif):
#         raise ValueError("motif length not equal to profile length")
#
#     n = len(motif)
#     pr = 1
#     for i in range(n):
#         pr *= profile.loc[motif[i], i]
#
#     return pr


# def ProfileMostProbableKmer(dna, profile):
#     n = len(dna)
#     k = profile.shape[1]
#     bestKmer = dna[0: k]
#     bestPr = Pr(bestKmer, profile)
#
#     for i in range(n - k + 1):
#         subdna = dna[i: i + k]
#         pr = Pr(subdna, profile)
#
#         if pr > bestPr:
#             bestPr = pr
#             bestKmer = subdna
#
#     return bestKmer


def Count(dnas, pseudo=False):
    n = len(dnas[0])
    if pseudo:
        count = {'A': [1]*n, 'C': [1]*n, 'G': [1]*n, 'T': [1]*n}
    else:
        count = {'A': [0]*n, 'C': [0]*n, 'G': [0]*n, 'T': [0]*n}

    for dna in dnas:
        for i in range(n):
            count[dna[i]][i] += 1

    return count

def Profile(dnas, pseudo=False):
    count = Count(dnas, pseudo)
    n = len(dnas[0])
    t = sum([count[nucleo][0] for nucleo in "ACGT"])

    for nucleo in "AGCT":
        for i in range(n):
            count[nucleo][i] /= t

    return count

def Pr(dna, profile):
    n = len(dna)
    assert n == len(profile['A']), "Profile length should be equal to DNA length"

    pr = 1
    for i, nucleo in enumerate(dna):
        pr *= profile[nucleo][i]

    return pr

def ProfileMostProbableKmer(dna, profile):
    k = len(profile['A'])
    n = len(dna)

    bestKmer = dna[0: k]
    bestPr = Pr(bestKmer, profile)
    for i in range(n - k + 1):
        kmer = dna[i: i+k]
        pr = Pr(kmer, profile)
        if pr > bestPr:
            bestPr = pr
            bestKmer = kmer

    return bestKmer

def Score(dnas):
    n = len(dnas[0])
    count = Count(dnas)

    score = 0
    for i in range(n):
        freqs = [count[nucleo][i] for nucleo in "ACGT"]
        maxFreq = max(freqs)
        score += sum(filter(lambda x: x != maxFreq, freqs))

    return score

def GreedyMotifSearch(dnas, k, pseudo=False):
    n = len(dnas[0])
    t = len(dnas)
    bestMotifs = [dna[0: k] for dna in dnas]

    for i in range(n - k + 1):
        print(f"Kmer {i} of {n - k + 1}")
        motifs = [dnas[0][i: i + k]]
        profile = Profile(motifs, pseudo=pseudo)

        for dna in dnas[1:]:
            motifs.append(ProfileMostProbableKmer(dna, profile))
            profile = Profile(motifs, pseudo=pseudo)
        if Score(motifs) < Score(bestMotifs):
            bestMotifs = motifs

    return bestMotifs

if __name__ == "__main__":
    pass
    # with open("week 3 datasets/dataset_160_9.txt") as file:
    #     k, t = map(int, file.readline().strip().split(" "))
    #     print(k)
    #     dnas = file.readline().strip().split(" ")
    #     print(dnas)
    #
    #     print(*GreedyMotifSearch(dnas, k, pseudo=True))

    # print(GreedyMotifSearch("GGCGTTCAGGCA AAGAATCAGTCA CAAGGAGTTCGC CACGTCAATCAC CAATAATATTCG".split(" "), 3,
    #                         pseudo=True))
    #
    # with open("week 3 datasets/dataset_5164_1.txt") as file:
    #     pattern = file.readline().strip()
    #     dnas = file.readline().strip().split(" ")
    #
    #     print(sum([distance(pattern, dna) for dna in dnas]))
    # print()


