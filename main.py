import numpy as np
from matplotlib import pyplot as plt

# frequency table of substrings of length `k` within `text`
def FrequencyTable(text, k):
    n = len(text)
    freqMap = {}

    for i in range(0, n - k + 1):
        pattern = text[i:i + k]
        if pattern in freqMap:
            freqMap[pattern] += 1
        else:
            freqMap[pattern] = 1

    return freqMap

# returns most frequent k-mer in `text`
def FrequentWords(text, k):
    n = len(text)
    freqPatterns = set()

    freqMap = FrequencyTable(text, k)

    maxCount = max(freqMap.values())
    for k, v in freqMap.items():
        if v == maxCount:
            freqPatterns.add(k)

    return freqPatterns

# returns amount `pattern` is in `text`
def PatternCount(text, pattern):
    n = len(text)
    k = len(pattern)
    count = 0
    for i in range(0, n - k + 1):
        substr = text[i:i+k]
        if substr == pattern:
            count += 1

    return count

def ReverseComplement(text):
    # Complement
    transTable = text.maketrans("atgcATGC", "tacgTACG")
    text = text.translate(transTable)
    return text[::-1]

# returns indices where `pattern` is in within `genome`
def PatternMatch(pattern, genome):
    n = len(genome)
    k = len(pattern)
    indices = []
    for i in range(0, n-k+1):
        subgenome = genome[i: i+k]
        if subgenome == pattern:
            indices.append(i)

    return indices

# find (L, t) clumps of length k in text
def FindClumps(text, k, L, t):
    patterns = set()
    n = len(text)
    for i in range(0, n-L+1):
        window = text[i:i+L]
        freqMap = FrequencyTable(window, k)

        for key, val in freqMap.items():
            if val >= t:
                patterns.add(key)

        if i%10000 == 0:
            print(f"{i} of {n-L+1}")
            print(*patterns, sep=" ")

    return patterns

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    with open("E_coli.txt") as lines_:
        genome = lines_.readline()
        print(*MinSkew(genome), sep=" ")





