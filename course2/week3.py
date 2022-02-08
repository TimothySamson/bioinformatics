amino = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"XXX", "UAG":"XXX",
    "UGU":"C", "UGC":"C", "UGA":"XXX", "UGG":"W",
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

def DnaToRna(dna):
    transTable = dna.maketrans("tTuU", "uUuU")
    text = dna.translate(transTable)
    return text

def ReverseComplementRna(rna):
    transTable = rna.maketrans("augcAUGC", "uacgUACG")
    text = rna.translate(transTable)
    return text[::-1]

def ProteinTranslate(string, stopAsXXX=False):
  codons = []
  for i in range(0, len(string), 3):
    codon = amino[string[i: i+3]]
    if codon == "XXX" and not stopAsXXX:
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

if __name__ == "__main__":
  # with open("course2/dataset_96_4.txt") as file:
  #   string = file.readline().strip()
  #   print(ProteinTranslate(string))

  # with open("dataset_96_7.txt") as file:
  #     dna = file.readline().strip()
  #     peptide = file.readline().strip()
  #     print(*PeptideEncoding(dna, peptide), sep="\n")

  with open("Bacillus_brevis.txt") as file:
      dna = "".join([line.strip() for line in file.readlines()])
      print(PeptideEncoding(dna, "VKLFPWFNQY"))
