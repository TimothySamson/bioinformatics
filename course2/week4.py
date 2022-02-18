import numpy as np
# import pandas as pd
from week3 import Spectrum, aminoMass, PeptideMass

def PeptideScore(peptide, spectrum, cyclic=True):
    theoretical = np.array(list(Spectrum(peptide, cyclic=cyclic)))
    a, b = np.unique(theoretical, return_counts=True)
    c, d = np.unique(spectrum, return_counts=True)
    theoreticalCounts = dict(zip(a, b))
    experimentalCounts = dict(zip(c, d))

    return sum([min(freq, theoreticalCounts[weight]) for weight, freq in experimentalCounts.items()
                                                     if weight in theoreticalCounts])

def Trim(leaderboard, spectrum, N):
    leaderScores = np.array([(peptide, PeptideScore(peptide, spectrum, cyclic=False)) for peptide in leaderboard],
                            dtype=[("peptide", object), ("score", int)])
    return [x[0] for x in np.sort(leaderScores, order="score")[::-1][:N]]

def LeaderboardCyclopeptideSequencing(spectrum, N):
    leaderboard = [""]
    finalists = []

    realMass = spectrum[-1]
    aminos = list(aminoMass.keys())

    while leaderboard:
        leaderboard = [
            previous + acid for acid in aminos for previous in leaderboard
        ]

        leaderboard = Trim(leaderboard, spectrum, N)
        print(f"{max([PeptideMass(peptide) for peptide in leaderboard])} of {realMass}")

        for peptide in leaderboard.copy():
            if PeptideMass(peptide) > realMass:
                leaderboard.remove(peptide)
            elif PeptideMass(peptide) == realMass:
                finalists.append(peptide)
                leaderboard.remove(peptide)

    return finalists



if __name__ == "__main__":
    # with open("dataset_4913_1.txt") as file:
    #     peptide = file.readline().strip()
    #     spectrum = np.array(file.readline().strip(" ").split(" ")).astype(int)
    #     print(PeptideScore(peptide, spectrum, cyclic=False))

    # with open("dataset_4913_3.txt") as file:
    #     leaderboard = file.readline().strip().split(" ")
    #     spectrum = np.array(file.readline().strip().split(" ")).astype(int)
    #     N = int(file.readline().strip())
    #     print(*[x[0] for x in Trim(leaderboard, spectrum, N)])

    with open("dataset_102_8.txt") as file:
        N = int(file.readline())
        spectrum = np.array(file.readline().strip().split(" ")).astype(int)

        res = LeaderboardCyclopeptideSequencing(spectrum, N)
        print(res)
        weights = set([(aminoMass[acid] for acid in peptide) for peptide in res])
        for weight in weights:
            print(list(weight))
            print(*list(weight), end=" ", sep="-")


