import copy
import sys
import numpy as np
import simplegraphs as sg

MAX_WIDTH = 2000
MAX_HEIGHT = 2000


def ifRounding(matrix):
    for row in matrix:
        if sum(row) % 10 != 0:
            return False
    Transpose = np.transpose(matrix)
    for row in Transpose:
        if sum(row) % 10 != 0:
            return False
    return True


def roundedMatrix(Gf, matrix, ans):
    T = np.transpose(matrix)
    for c, rt in enumerate(T):
        rt = str(rt)
        if rt not in Gf['adj']['t']:
            continue
        for r, l in enumerate(matrix):
            l = str(l)
            if l in Gf['adj'][rt]:
                ans[r][c] = ans[r][c] + (10 - (matrix[r][c] % 10))
    for r in range(len(ans)):
        for c in range(len(ans[0])):
            if ans[r][c] % 10 != 0:
                ans[r][c] = ans[r][c] - ans[r][c] % 10
    return ans


def augmentingPath(G, s, t, method="EdmonsonKarp2"):
    if method == "FordFulkerson":
        parent = sg.DFS(G)[2]
        if parent[t] is None:
            return -1, G
    elif method == "EdmonsonKarp2":
        parent = sg.BFS(G, s)[1]
        if t not in parent:
            return -1, G
    temp = t
    minimumCut = 1e15
    while temp != s:
        tempParent = parent[temp]
        capacity = G['adj'][tempParent][temp]
        if capacity < minimumCut:
            minimumCut = capacity
        temp = tempParent
    temp = t
    while temp != s:
        tempParent = parent[temp]
        G['adj'][tempParent][temp] = G['adj'][tempParent][temp] - minimumCut
        if G['adj'][tempParent][temp] == 0:
            sg.delEdge(G, tempParent, temp)
        if temp in G['adj']:
            if tempParent not in G['adj'][temp]:
                sg.addDirEdge(G, temp, tempParent, label=minimumCut)
        else:
            G['adj'][temp][tempParent] = G['adj'][temp][tempParent] + minimumCut
        temp = tempParent
    return minimumCut, G


def newGoldGraph(coords):
    l, r = set(), set()
    G = {'n': 0, 'm': 0, 'adj': {}}
    for array in coords:
        if array[0] % 2 == 1 and array[1] % 2 == 1:
            l.add(f'{array[0]} {array[1]}')
        elif array[0] % 2 == 0 and array[1] % 2 == 0:
            l.add(f'{array[0]} {array[1]}')
        else:
            r.add(f'{array[0]} {array[1]}')
    if len(r) != len(l):
        return False, G
    def addEdge(Gf, n, m):
        sg.addDirEdge(Gf, 's', n, 1)
        sg.addDirEdge(Gf, n, m, 1)
        sg.addDirEdge(Gf, m, 't', 1)
        return Gf
    for x in l:
        node1, node2 = x.split()
        node1, node2 = int(node1), int(node2)
        if f'{str(node1 + 1)} {str(node2)}' in r:
            G = addEdge(G, f'{node1} {node2}', f'{node1 + 1} {node2}')
        if f'{str(node1)} {str(node2 + 1)}' in r:
            G = addEdge(G, f'{node1} {node2}', f'{node1} {node2 + 1}')
        if f'{str(node1 - 1)} {str(node2)}' in r:
            G = addEdge(G, f'{node1} {node2}', f'{node1 - 1} {node2}')
        if f'{str(node1)} {str(node2 - 1)}' in r:
            G = addEdge(G, f'{node1} {node2}', f'{node1} {node2 - 1}')
    return True, G


def newRoundingGraph(matrix):
    G = {'n': 0, 'm': 0, 'adj': {}}
    for r in matrix:
        capacity = 0
        for c in r:
            capacity = capacity + c % 10
        r = str(r)
        sg.addDirEdge(G, 's', r, capacity)
    T = np.transpose(matrix)
    for r in T:
        capacity = 0
        for c in r:
            capacity = capacity + c % 10
        r = str(r)
        sg.addDirEdge(G, r, 't', capacity)
    for r in range(len(matrix)):
        for c in range(len(matrix[r])):
            if matrix[r][c] % 10 != 0:
                edge = str(T[c])
                sg.addDirEdge(G, str(matrix[r]), edge, 10)
    return G


def gold(coords):
    ans = []
    if len(coords) % 2 != 0:
        print('No solution exists')
        return 'No solution exists'
    present, G = newGoldGraph(coords)
    if not present:
        print('No solution exists')
        return 'No solution exists'
    Gf = maxflow(G, 's', 't', method="EdmonsonKarp2")
    if len(coords) / 2 > len(Gf['adj']['t']):
        print('No solution exists')
        return 'No solution exists'
    for r in Gf['adj']['t']:
        l = list(Gf['adj'][r].keys())[0]
        ans.append(f'{l} --> {r}')
    for i in ans:
        print(i)
    if len(ans) == 0:
        print('No solution exists')
        return 'No solution exists'
    else:
        return


def rounding(matrix):
    if not ifRounding(matrix):
        print('No solution exists')
        return matrix
    Gf = maxflow(newRoundingGraph(matrix), 's', 't', method='FordFulkerson')
    ans = roundedMatrix(Gf, matrix, copy.deepcopy(matrix))
    return ans


def maxflow(G, s, t, method):
    Gf = sg.copyGraph(G)
    maxFlow = 0
    while 1:
        minimumCut, Gf = augmentingPath(Gf, s, t, method)
        if minimumCut != -1:
            maxFlow = maxFlow + minimumCut
        else:
            break
    return Gf


def main(args=[]):
    if len(args) < 2:
        print("Too few arguments! There should be at least 4.")
        print("flows.py <cmd> <file>")
        return
    task = args[0]
    if task == "gold":
        coords = read_input(args[1])
        gold(coords)
    elif task == "rounding":
        matrix = read_input(args[1])
        nm = rounding(matrix)
        if compare_matrix(matrix, nm):
            print_matrix(nm)
    elif task == "maxflow":
        graph_file = args[1]
        s = int(args[2])
        t = int(args[3])
        G = sg.readGraph(graph_file)
        flow = maxflow(G, s, t)
        print(flow)
    return


def read_input(filename):
    with open(filename, 'r') as f:
        blocks = [[int(x) for x in s.split()] for s in f.read().splitlines()]
    return blocks


def print_matrix(matrix):
    for r in matrix:
        print(*r)


def compare_matrix(m1, m2):
    r1 = len(m1)
    r2 = len(m2)
    c1 = len(m1[0])
    c2 = len(m2[0])
    if r1 != r2 or c1 != c2:
        print('Sizes are different')
        return False
    for ri in range(0, r1):
        rs1 = sum(m1[ri])
        rs2 = sum(m2[ri])
        if rs1 != rs2:
            print('Row sum {ri} differ: {rs1} != {rs2}')
            return False
    for cj in range(0, c1):
        cs1 = 0
        cs2 = 0
        for i in range(0, r1):
            cs1 += m1[i][cj]
            cs2 += m2[i][cj]
        if cs1 != cs2:
            print('Col sum {cj} differ: {cs1} != {cs2}')
            return False
    return True


if __name__ == "__main__":
    main(sys.argv[1:])