# returns the overlap adjacency list (hamiltonian)
def OverlapAdjacencyList(dnas):
    adjacencyList = {}
    for dnaKey in dnas:
        if dnaKey not in adjacencyList:
            adjacencyList[dnaKey] = []
        for dna in dnas:
            if dna == dnaKey:
                continue
            if dnaKey[1:] == dna[:-1]:
                adjacencyList[dnaKey].append(dna)

    return adjacencyList


def PrintAdjacencyList(adjList):
    for key, value in adjList.items():
        if value:
            print(f"{key} -> ", end="")
            print(*value, sep=",")


def HamiltonianPaths(adjList, start):
    paths = []
    def HamiltonianPathInner(start, path=()):
        assert start in adjList, "start node not in adjacency list"

        path = path + (start, )
        print(path)
        if len(path) == len(adjList.keys()):
            paths.append(path)


        for node in adjList[start]:
            if node not in path:
                HamiltonianPathInner(node, path)

    HamiltonianPathInner(start)
    return paths


def nbits(n):
    if n == 1:
        return ["0", "1"]

    prevBits = nbits(n - 1)
    return ["0" + bit for bit in prevBits] + ["1" + bit for bit in prevBits]

# Returns all k-mers in text
def Composition(text, k):
    return [text[i: i + k] for i in range(len(text) - k + 1)]

# Make DeBruijn adjacency list from kmers
def KmerToDeBruijn(kmers):
    graph = {kmer[:-1]: [] for kmer in kmers}
    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        graph[prefix].append(suffix)
        if suffix not in graph:
            graph[suffix] = []


    return graph

# Make debruijn graph from a text
def DeBruijnGraph(text, k):
    n = len(text)
    graph = {}
    for i in range(n - k + 2):
        key = text[i: i+k-1]
        if key in graph:
            continue

        graph[key] = []
        for j in range(n-k+1):
            word = text[j: j+k-1]
            if word == key:
                graph[key].append(text[j+1: j+k])

    return graph

def PathToGenome(dnas):
    n = len(dnas[0])
    for i in range(n - 1):
        assert dnas[i][1:] == dnas[i + 1][:-1], "Not a path"

    return dnas[0] + "".join([dna[-1] for dna in dnas[1:]])

if __name__ == "__main__":
    # with open("course2/dataset_197_3 (1).txt") as file:
    #   k = int(file.readline().strip())
    #   text = file.readline().strip()
    #
    #   print(*Composition(text, k), sep="\n")
    #
    # with open("course2/dataset_198_3.txt") as file:
    #   dnas = [line.strip() for line in file.readlines()]
    #   print(dnas)
    #   print(PathToGenome(dnas))

    # with open("course2/dataset_198_10.txt") as file:
    #   dnas = [line.strip() for line in file.readlines()]
    #   PrintAdjacencyList(OverlapAdjacencyList(dnas))

    adjList = {
        0: [1],
        1: [2],
        2: [3, 1],
        3: [0, 1]
    }

    # print(HamiltonianPath(adjList, 0))
    # print(nbits(4))
    adj = OverlapAdjacencyList(nbits(4))
    PrintAdjacencyList(adj)

    print(*[PathToGenome(path) for path in HamiltonianPaths(adj, "0000")], sep="\n")