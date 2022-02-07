from course2.week1 import PrintAdjacencyList, KmerToDeBruijn, PathToGenome, nbits
import re
from copy import deepcopy
import itertools

def ParseToAdjList(pairs):
    adjList = {}
    for pair in pairs:
        match = re.search("(.+) -> (.+)", pair)
        key = match.group(1)
        value = match.group(2).split(",")
        adjList[key] = value

    return adjList

# Randomly walks along graph and tries to make a cycle. Modifies its argument
def EulerianCycleInner(adjList, start):
    path = []
    curNode = start
    while adjList[curNode]:
        nextNodes = adjList[curNode]
        curNode = nextNodes.pop()
        path.append(curNode)
    return path


def EulerianCycle(adjList, start):
    adjList = deepcopy(adjList)
    path = EulerianCycleInner(adjList, start)
    while True:
        toBreak = True
        for i, node in enumerate(path):
            if adjList[node]:
                path[i+1: i+1] = EulerianCycleInner(adjList, node)
                toBreak = False

        if toBreak:
            break

    return [start] + path

def isEulerCycle(adjList, path):
    adjList = deepcopy(adjList)
    curNode = path[0]

    for i, nextNode in enumerate(path[1:]):
        assert nextNode in adjList[curNode], f"{nextNode} is not a valid edge from {curNode} in place {i+1}"
        adjList[curNode].remove(nextNode)
        curNode = nextNode

    return not any(adjList[node] for node in path)

# for each node, a dict with key "out" & "in"
def DegreeList(adjList):
    degree = {}

    for node, neighbors in adjList.items():
        # initialize degree dict
        for vertex in neighbors + [node]:
            if vertex not in degree:
                degree[vertex] = {"in": 0, "out": 0}

        degree[node]["out"] = len(neighbors)
        for neighbor in neighbors:
            degree[neighbor]["in"] += 1

    return degree

def EulerianPath(adjList):
    # Figure out end and start node
    adjList = deepcopy(adjList)
    degree = DegreeList(adjList)
    for node in adjList.keys():
        diff = degree[node]["out"] - degree[node]["in"]
        # print(node, diff)
        if diff == 1:
            start = node
        if diff == -1:
            end = node

    assert "start" in locals(), "no start node"
    assert "end" in locals(), "no end node"

    # add arbitrary edge
    adjList[end].append(start)

    # rearrange cycle
    cycle = EulerianCycle(adjList, start)
    startIndices = [i for i, node in enumerate(cycle) if node == start]
    cycle.pop(0)

    loops = [cycle[startIndices[i]: startIndices[i+1]] for i in range(len(startIndices) - 1)]

    lastLoopIndex = next(i for i, loop in enumerate(loops) if loop[-2:] == [end, start])
    loops[lastLoopIndex], loops[-1] = loops[-1], loops[lastLoopIndex]

    cycle = [start] + sum(loops, [])
    cycle.pop()
    return cycle

def CircularUniversalString(k):
    bits = nbits(k)
    adjList = KmerToDeBruijn(bits)
    return PathToGenome(EulerianCycle(adjList, "0"*(k-1)))[:-k+1]

def PairedPathToGenome(path, d):
    k = len(path[0][0])
    leftReads = [read[0] for read in path]
    rightReads = [read[1] for read in path]
    prefix = leftReads[0] + "".join([read[-1] for read in leftReads[1:]])
    suffix = rightReads[0] + "".join([read[-1] for read in rightReads[1:]])

    leftOverlap = prefix[k+d:]
    rightOverlap = suffix[:-k-d]

    print(prefix)
    print(" " * (k + d) + suffix)

    assert leftOverlap == rightOverlap, "no genome spelled by gapped pairs"

    return prefix + suffix[len(leftOverlap):]

def ReadPairAdjacencyList(reads):
    adjList = {(read[0][:-1], read[1][:-1]): [] for read in reads}

    for read in reads:
        suffix = (read[0][1:], read[1][1:])
        prefix = (read[0][:-1], read[1][:-1])
        adjList[prefix].append(suffix)
        if suffix not in adjList:
            adjList[suffix] = []

    return adjList

# def MaximalNonBranchingPaths(graph):
#     graph = deepcopy(graph)
#     degree = DegreeList(graph)
#
#     for node in graph.keys():
#         if degree[node] != {"out": 1, "in", 1}:
#
#



if __name__ == "__main__":
    # with open("dataset_203_6.txt") as file:
    #     adjList = ParseToAdjList([line.strip() for line in file.readlines()])
    #     res = EulerianPath(adjList)
    #
    #     print(*res, sep="->")

    # with open("dataset_203_7.txt") as file:
    #     adjList = KmerToDeBruijn([line.strip() for line in file.readlines()])
    #     res = EulerianPath(adjList)
    #     print(PathToGenome(res))

    # with open("dataset_6206_4.txt") as file:
    #     path = [read.strip().split("|") for read in file.readlines()]
    #     print(PairedPathToGenome(path, 200))

    # with open("dataset_204_16.txt") as file:
    #     pairs = [read.strip().split("|") for read in file.readlines()]
    #     d = 200
    #     adjList = ReadPairAdjacencyList(pairs, d)
    #     print(PairedPathToGenome(EulerianPath(adjList), d+1))

    with open("input.txt") as file:
        pairs = [read.strip().split("|") for read in file.readlines()]
        adjList = ReadPairAdjacencyList(pairs)
        PrintAdjacencyList(adjList)
        print(PairedPathToGenome(EulerianPath(adjList), 2))







