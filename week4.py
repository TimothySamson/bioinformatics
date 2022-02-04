import pandas as pd
import numpy as np
import random

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

def GibbsSamplerIter(dnas, k, N=1000):
    n = len(dnas[0])
    t = len(dnas)
    randIndices = np.random.randint(n - k + 1, size=t)

    bestMotifs = [dna[randIndices[i]: randIndices[i] + k] for (i, dna) in enumerate(dnas)]
    bestScore = Score(bestMotifs)

    motifs = bestMotifs

    for i in range(N):
        # Removes random motif
        popIndex = np.random.randint(t)
        motifs.pop(popIndex)

        # Makes profile of the rest of motifs, and makes the probability distribution
        removedIndexProfile = Profile(motifs, pseudo=True)
        dna = dnas[popIndex]
        prs = [Pr(dna[i: i+k], removedIndexProfile) for i in range(n-k+1)]
        normPrs = list(map(lambda x: x / sum(prs), prs))

        # Choose from a probability distribution, and insert it to motifs
        motifIndex = np.random.choice(n-k+1, p=normPrs)
        newMotif = dna[motifIndex: motifIndex + k]
        motifs.insert(popIndex, newMotif)

        score = Score(motifs)

        if score < bestScore:
            # print(f"{score}: new best motifs: ", *motifs)
            bestMotifs = motifs
            bestScore = score

    return bestMotifs

def GibbsSampler(dnas, k, N=2000, starts=20):
    bestMotifs = GibbsSamplerIter(dnas, k, N)
    bestScore = Score(bestMotifs)

    for i in range(starts):
        motifs = GibbsSamplerIter(dnas, k, N)
        score = Score(motifs)

        print(f"{score}: gibbs {i} of {starts} (best: {bestScore})")

        if score < bestScore:
            bestMotifs = motifs
            bestScore = score

    return bestMotifs


def RandomizedMotifSearchIter(dnas, k):
    n = len(dnas[0])
    t = len(dnas)

    randIndices = np.random.randint(n - k + 1, size=t)
    bestMotifs = [dna[randIndices[i]: randIndices[i] + k] for (i, dna) in enumerate(dnas)]

    while True:
        profile = Profile(bestMotifs, pseudo=True)
        randMotifs = [ProfileMostProbableKmer(dna, profile) for dna in dnas]

        if Score(bestMotifs) <= Score(randMotifs):
            return bestMotifs

        bestMotifs = randMotifs


def RandomizedMotifSearch(dnas, k, n=1000):
    bestMotifs = RandomizedMotifSearchIter(dnas, k)
    print(len(dnas[0]))

    for i in range(n):
        randMotifs = RandomizedMotifSearchIter(dnas, k)
        print(f"{Score(randMotifs)}: Random Iter {i} of {n}")
        if Score(randMotifs) < Score(bestMotifs):
            bestMotifs = randMotifs

    return bestMotifs


if __name__ == "__main__":
    # print(RandomizedMotifSearch("CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG TAGTACCGAGACCGAAAGAAGTATACAGGCGT TAGATCAAGTTTCAGGTGCACGTCGGTGAACC AATCCACCAGCTCCACGTGCAATGTTGGCCTA".split(" "), 8))

    # with open("dataset_161_5.txt") as file:
    #     k, t = map(int, file.readline().strip().split(" "))
    #     dnas = file.readline().strip().split(" ")
    #     A = RandomizedMotifSearch(dnas, k)
    #     print(*A)
    #     print(Score(A))

    # profile = {
    #     'A': [0.2, 0.2, 0.3, 0.2, 0.3],
    #     'C': [0.4, 0.3, 0.1, 0.5, 0.1],
    #     'G': [0.3, 0.3, 0.5, 0.2, 0.4],
    #     'T': [0.1, 0.2, 0.1, 0.1, 0.2]
    # }

    # print(Score("TCTCGGGG CCAAGGTG TACAGGCG TTCAGGTG TCCACGTG".split(" ")))
    # print(Count("TCTCGGGG CCAAGGTG TACAGGCG TTCAGGTG TCCACGTG".split(" ")))
    # print(Profile("TCTCGGGG CCAAGGTG TACAGGCG TTCAGGTG TCCACGTG".split(" "), pseudo=True))

    # print(ProfileMostProbableKmer("ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT", profile))

    # with open("dataset_163_4.txt") as file:
    #     k, t, N = map(int, file.readline().strip().split(" "))
    #     dnas = file.readline().strip().split(" ")
    #
    #     print(*GibbsSampler(dnas, k, N=2000, starts=20))


    # dnas = "CGCCCCTCTCGGGGGTGTTCAGTAACCGGCCA GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG TAGTACCGAGACCGAAAGAAGTATACAGGCGT TAGATCAAGTTTCAGGTGCACGTCGGTGAACC AATCCACCAGCTCCACGTGCAATGTTGGCCTA".split(" ")
    #
    # res = GibbsSampler(dnas, 8, 1000)
    # print(*res)
    # ans = "TCTCGGGG CCAAGGTG TACAGGCG TTCAGGTG TCCACGTG".split(" ")
    #
    # print("Score of ans: ", Score(ans))
    # print("Score of res: ", Score(res))
    #
    #

    dnas = "ATGAGGTC GCCCTAGA AAATAGAT TTGTGCTA".split(" ")
    motifs = "GTC CCC ATA GCT".split(" ")
    profile = Profile(motifs)
    nextMotifs = [ProfileMostProbableKmer(dna, profile) for dna in dnas]
    print(*nextMotifs)
