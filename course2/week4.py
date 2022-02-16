import numpy as np
# import pandas as pd
from week3 import Spectrum

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
    return np.sort(leaderScores, order="score")[::-1][:N]




if __name__ == "__main__":
    # with open("dataset_4913_1.txt") as file:
    #     peptide = file.readline().strip()
    #     spectrum = np.array(file.readline().strip(" ").split(" ")).astype(int)
    #     print(PeptideScore(peptide, spectrum, cyclic=False))

    with open("dataset_4913_3.txt") as file:
        leaderboard = file.readline().strip().split(" ")
        spectrum = np.array(file.readline().strip().split(" ")).astype(int)
        N = int(file.readline().strip())
        print(*[x[0] for x in Trim(leaderboard, spectrum, N)])

