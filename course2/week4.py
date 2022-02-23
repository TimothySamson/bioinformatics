import numpy as np
# import pandas as pd
from week3 import Spectrum, aminoMass

def PeptideScore(peptide, spectrum, cyclic=True):
    a, b = np.unique(list(Spectrum(peptide, cyclic=cyclic)),
                     return_counts=True)
    c, d = np.unique(spectrum, return_counts=True)
    theoreticalCounts = dict(zip(a, b))
    experimentalCounts = dict(zip(c, d))
    return sum([
        min(freq, theoreticalCounts[weight])
        for weight, freq in experimentalCounts.items()
        if weight in theoreticalCounts
    ])

def Trim(leaderboard, spectrum, N):
    if len(leaderboard) <= N:
        return leaderboard

    leaderboard = np.array(leaderboard)
    leaderScores = np.array(
        [PeptideScore(peptide, spectrum, cyclic=False) for peptide in leaderboard]
    )

    upToScore = leaderScores[np.argsort(-leaderScores)][N-1]
    return list(tuple(x) for x in leaderboard[leaderScores >= upToScore])


def LeaderboardCyclopeptideSequencing(spectrum, N, extended=False):
    leaderboard = [()]
    if not extended:
        aminos = list(set(aminoMass.values()))
    else:
        aminos = list(range(57, 200 + 1))

    finalists = []
    realMass = spectrum[-1]

    while leaderboard:
        leaderboard = [previous + (acid,) for acid in aminos for previous in leaderboard]

        leaderboard = Trim(leaderboard, spectrum, N)

        for peptide in leaderboard.copy():
            if sum(peptide) > realMass:
                leaderboard.remove(peptide)
            elif sum(peptide) == realMass:
                finalists.append(peptide)
        if leaderboard:
            print(f"{max([sum(peptide) for peptide in leaderboard])} of {realMass}")

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

    # with open("dataset_102_8.txt") as file:
    #     N = int(file.readline())
    #     spectrum = np.array(file.readline().strip().split(" ")).astype(int)
    #
    #     res = np.array(LeaderboardCyclopeptideSequencing(spectrum, N))
    #     scores = np.array([PeptideScore(peptide, spectrum, cyclic=True) for peptide in res])
    #     res = res[scores == max(scores)]
    #     for weight in res:
    #         print(*weight, end=" ", sep="-")

    # with open("dataset_103_2.txt") as file:
    #     N = int(file.readline().strip())
    #
    #     spectrum = np.array(file.readline().strip().split(" ")).astype(int)
    #
    #     for peptide in LeaderboardCyclopeptideSequencing(spectrum, N, extended=True):
    #         print(*peptide, sep="-", end=" ")

    # with open("input.txt") as file, open("dataset_103_2.txt") as spectrum:
    #     N = int(spectrum.readline().strip())
    #     spectrum = np.array(spectrum.readline().strip().split(" ")).astype(int)
    #
    #     peptides = np.array([tuple(map(int, peptide.split("-"))) for peptide in file.readline().strip().split(" ")])
    #     pepScores = [PeptideScore(peptide, spectrum, cyclic=False, extended=True) for peptide in peptides]
    #     pepScores = np.array(pepScores)
    #
    #     print(len(peptides[pepScores == max(pepScores)]))
    #     print(pepScores)

    spectrum = "0 97 99 113 114 115 128 128 147 147 163 186 227 241 242 244 244 256 260 261 262 283 291 309 330 333 340 347 385 388 389 390 390 405 435 447 485 487 503 504 518 544 552 575 577 584 599 608 631 632 650 651 653 672 690 691 717 738 745 770 779 804 818 819 827 835 837 875 892 892 917 932 932 933 934 965 982 989 1039 1060 1062 1078 1080 1081 1095 1136 1159 1175 1175 1194 1194 1208 1209 1223 1322"
    spectrum = np.array(spectrum.split(" ")).astype(int)

    peptides = np.array(LeaderboardCyclopeptideSequencing(spectrum, 1000))
    peptideScores = np.array([PeptideScore(peptide, spectrum) for peptide in peptides])

    peptides = peptides[peptideScores == max(peptideScores)]
    print(peptides)
    print(len(peptides))





