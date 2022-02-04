from course2.week1 import PrintAdjacencyList, KmerToDeBruijn, PathToGenome
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
    degree = DegreeList(adjList)
    for node in adjList.keys():
        diff = degree[node]["out"] - degree[node]["in"]
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


if __name__ == "__main__":
    # with open("dataset_203_6.txt") as file:
    #     adjList = ParseToAdjList([line.strip() for line in file.readlines()])
    #     res = EulerianPath(adjList)
    #
    #     print(*res, sep="->")

    with open("dataset_203_7.txt") as file:
        adjList = KmerToDeBruijn([line.strip() for line in file.readlines()])
        res = EulerianPath(adjList)
        print(PathToGenome(res))






