import itertools
from copy import deepcopy

amino = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"X", "UAG":"X",
    "UGU":"C", "UGC":"C", "UGA":"X", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"}

aminos = list(set(amino.values()))
aminoMass = {"G": 57, "A": 71 , "S": 87 , "P": 97 , "V": 99 , "T": 101, "C": 103, "I": 113, "L": 113, "N": 114,
        "D": 115,  "K": 128,  "Q": 128,  "E": 129,  "M": 131,  "H": 137,  "F": 147,  "R": 156,  "Y": 163,  "W": 186}





def DnaToRna(dna):
    transTable = dna.maketrans("tTuU", "uUuU")
    text = dna.translate(transTable)
    return text

def ReverseComplementRna(rna):
    transTable = rna.maketrans("augcAUGC", "uacgUACG")
    text = rna.translate(transTable)
    return text[::-1]

# Takes in RNA, returns protein
def ProteinTranslate(string, stopAsX=False):
  codons = []
  for i in range(0, len(string), 3):
    codon = amino[string[i: i+3]]
    if codon == "X" and not stopAsX:
      break
    codons.append(codon)
  return "".join(codons)

def PeptideEncoding(dna, peptide):
    subdnas = []
    n = len(dna)
    k = len(peptide)

    for i in range(n - k*3):
        # if i % 100000 == 0:
        #     print(f"{i} of {n-k*3}")
        kmer = dna[i:i + 3*k]
        trans = DnaToRna(kmer)
        if ProteinTranslate(trans) == peptide or ProteinTranslate(ReverseComplementRna(trans)) == peptide:
            print(f"found {kmer} at position {i}")
            subdnas.append(kmer)

    return subdnas

def CountPeptides(mass):
    aminoAcids = list(aminoMass.keys()).copy()
    aminoAcids.remove("L")
    aminoAcids.remove("Q")

    combs = {aminoMass[acid]: 1 for acid in aminoAcids}

    while any(map(lambda x: x < mass, combs.keys())):
        newCombs = {}
        print(combs)
        for partialMass, comb in combs.items():
            if partialMass >= mass:
                newCombs[partialMass] = comb
                continue
            for acid in aminoAcids:
                newMass = partialMass + aminoMass[acid]
                if newMass not in newCombs:
                    newCombs[newMass] = comb
                else:
                    newCombs[newMass] += comb
        combs = newCombs


    print(combs)
    return combs[mass]


def Spectrum(peptide, cyclic=True):
    n = len(peptide)

    yield 0

    # accessible until n
    prefixMass = [0] + list(itertools.accumulate([aminoMass[acid] for acid in peptide]))

    peptideMass = prefixMass[-1]
    for i in range(n):
        for j in range(i+1, n+1):
            middleMass = prefixMass[j] - prefixMass[i]

            yield middleMass

            if i != 0 and j != n and cyclic:
                yield peptideMass - middleMass

def isSubset(subset, parent):
  parent = deepcopy(parent)
  for elem in subset:
    if elem in parent:
      parent.remove(elem)
    else:
      return False
  return True

def PeptideMass(peptide):
  return sum([aminoMass[acid] for acid in peptide])

def CyclopeptideSequencing(spectrum):
    aminos = aminoMass.keys()
    aminos = [acid for acid in aminos if aminoMass[acid] in spectrum]

    candidatePeptides = [""]
    finalPeptides = []
    foundMatch = False
    finalMass = spectrum[-1]
    i = 0
    while not foundMatch:
        # Branch
        candidatePeptides = ([peptide + acid for acid in aminos
                              for peptide in candidatePeptides])

        i += 1
        print(f"{i}: {len(candidatePeptides)}")
        # Bound
        for peptide in candidatePeptides.copy():

            if PeptideMass(peptide) == finalMass:
                spectrumList = list(Spectrum(peptide))
                spectrumList.sort()
                if spectrumList == spectrum:
                    finalPeptides.append(peptide)
                    foundMatch = True

            elif not isSubset(Spectrum(peptide, cyclic=False), spectrum):
                candidatePeptides.remove(peptide)
    return finalPeptides


if __name__ == "__main__":
    # with open("course2/dataset_96_4.txt") as file:
    #   string = file.readline().strip()
    #   print(ProteinTranslate(string))

    # with open("dataset_96_7.txt") as file:
    #     dna = file.readline().strip()
    #     peptide = file.readline().strip()
    #     print(*PeptideEncoding(dna, peptide), sep="\n")

    # with open("Bacillus_brevis.txt") as file:
    #     dna = "".join([line.strip() for line in file.readlines()])
    #     print(PeptideEncoding(dna, "VKLFPWFNQY"))
    #   peptide = "KCPRFERQSWWDM"
    #   print(*CyclicSpectrum(peptide), sep=" ")


    # a = "0 97 99 113 114 115 128 129 137 147 147 186 212 " \
    #     "229 234 236 241 244 257 261 283 294 333 333 340 " \
    #     "349 358 370 372 376 408 420 430 446 469 477 480 " \
    #     "485 486 505 519 523 567 574 577 584 594 599 606 " \
    #     "632 633 652 666 691 698 703 709 714 721 746 760 " \
    #     "779 780 806 813 818 828 835 838 845 889 893 907 " \
    #     "926 927 932 935 943 966 982 992 1004 1036 1040 " \
    #     "1042 1054 1063 1072 1079 1079 1118 1129 1151 " \
    #     "1155 1168 1171 1176 1178 1183 1200 1226 1265 1265 " \
    #     "1275 1283 1284 1297 1298 1299 1313 1315 1412".split(" ")
    # with open("dataset_100_6.txt") as file:
    #     a = CyclopeptideSequencing(list(map(int, file.readline().strip().split(" "))))
    #     b = [(aminoMass[acid] for acid in peptide) for peptide in a]
    #     b = set(b)
    #
    #     for weightSet in b:
    #         print(*weightSet, end=" ", sep="-")
    # print([[PeptideMass(acid) for acid in peptide]
    #                           for peptide in CyclopeptideSequencing(list(map(int, a)))])

    # for a in Spectrum("NQEL"):
    #     print(a)

    # a = list(Spectrum("NQEL"))
    # a.sort()
    # print(a)

    # a = list(Spectrum("IQADNKIERIDRCTEQPVIHVFNKANRQAWMIIMGF", cyclic=False))
    # a.sort()
    # print(*a)

    print(ProteinTranslate("CCAAGAACAGAUAUCAAU"))

    a = list(map(int, "0 71 99 101 103 128 129 199 200 204 227 230 231 298 303 328 330 332 333".split(" ")))
    print([isSubset(Spectrum(peptide, cyclic=False), a) for peptide in [
        "TVQ", "QCV", "ETC", "TCQ", "TCE", "AVQ"
    ]])


