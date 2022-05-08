############################################################
# A simple graphs package
# Version 2.0
# October 2021
# Adam Smith
############################################################



import sys
import heapq
import queue
import numpy as np


def readGraph(input_file):
    '''This procedure takes the name of a text file describing a directed or undirected graph and returns a data structure in memory. 
    The file is structured as follows: each line lists an edge as a pair u,v, for an unweighted graph, or triples u,v,w for a weighted graph, where w represents the edge's weight as a float.
    The data structure it returns is a dictionary with keys "n", "m", "adj" (strings). The values for "n" and "m" are the number of nodes and edges in the graph (respectively). Edges are counted as in a directed graph (so each undirected edge counts twice). The value for "adj" is a dictionary-of-dictionaries adjacency structure. 
    '''
    with open(input_file, 'r') as f:
        raw = [[float(x) for x in s.split(',')] for s in f.read().splitlines()]
    G = emptyGraph(0)
    for edge in raw:
        if len(edge)==2:
            addDirEdge(G, int(edge[0]), int(edge[1]), label= True)
        else:
            addDirEdge(G, int(edge[0]), int(edge[1]), label= int(edge[2]))
    return G



def writeGraph(G, output_file):
    # G is a dictionary with keys "n", "m", "adj" representing a
    # weighted or unweighted graph
    with open(output_file, 'w') as f:
        for u in G["adj"].keys():
            for v in G["adj"][u]:
                if G["adj"][u][v] == True:
                    f.write("{}, {}\n".format(u,v))
                else:
                    f.write("{}, {}, {}\n".format(u,v, G["adj"][u][v] ))
    return



############################################################
# Functions for Basic Manipulations
############################################################

def copyGraph(G):
    # This will create a fresh copy of G in memory. Useful in case you
    # want to make changes but keep the original graph around.
    newG = {}
    newG["n"] = G["n"]
    newG["m"] = G["m"]
    newG["adj"] = {}
    for u in G["adj"]:
        newG["adj"][u] = {} # create a fresh dict for u's adjacency list
        for v in  G["adj"][u]:
            newG["adj"][u][v] = G["adj"][u][v] # copy whatever value was stored in G
    return newG

def reverseGraph(G):
    # This will create a fresh copy of G in memory. Useful in case you want to make changes but keep the original graph around.
    newG = {}
    newG["n"] = G["n"]
    newG["m"] = G["m"]
    newG["adj"] = {}
    for u in G["adj"]:
        newG["adj"][u] = {} # create a fresh dict for u's adjacency list
    for u in G["adj"]:
        for v in  G["adj"][u]:
            # copy whatever value was stored in G
            # but reverse the direction of the edge
            newG["adj"][v][u] = G["adj"][u][v] 
    return newG

def degree(G, u):
    return len(G["adj"][u])


def addNode(G,  x):
    # add a new node with name x
    # no effect if x is already in the graph
    if not(x in G["adj"]):
        G["adj"][x] = {} #create a new adjacency dict for x
        G["n"] = G["n"] + 1
    return
        
def addUndirEdge(G,  u, v, label = True):
    # add a new undirected edge from u to v with name x
    addNode(G, u)
    addNode(G, v)
    G["adj"][u][v] = label
    G["adj"][v][u] = label
    G["m"] = G["m"]  + 2
    return

def addDirEdge(G,  u, v, label = True):
    # add a new *directed* edge from u to v with name x
    addNode(G, u)
    addNode(G, v)
    G["adj"][u][v] = label
    G["m"] = G["m"]  + 1
    return

def delUndirEdge(G, u, v):
    # This will remove the edges (u,v) and (v,u) from G
    # It will throw an error if u,v or the edge do not exist.
    assert u in G["adj"] and v in G["adj"]
    assert v in G["adj"][u] and u in  G["adj"][v]
    del G["adj"][u][v]
    del G["adj"][v][u]
    G["m"] = G["m"] - 2
    return

def delEdge(G,u,v):
    assert u in G["adj"] and v in G["adj"]
    assert v in G["adj"][u]
    del G["adj"][u][v]
    G["m"] = G["m"] - 1
        
def equal(G, H):
    # determine if two graphs have the same set of directed edges
    if G["n"] != H["n"] or G["m"] != H["m"] or len(G["adj"]) != len(H["adj"]):
        return False
    for u in G["adj"]:
        if u not in H["adj"] or len(G["adj"][u]) != len(H["adj"][u]):
            return False
        for v in G["adj"][u]:
            if v not in H["adj"][u] or G["adj"][u][v] != H["adj"][u][v]:
                return False
    return True

def makeUndirected(G):
    # G is a dictionary with keys "n", "m", "adj" representing an unweighted directed graph
    # This function will modify G to  ensure that every edge (u,v) also appears as (v,u) in the adjacency structure.
    # If G is already undirected, this will have no effect.
    adj = G["adj"]
    for u in adj:
        for v in adj[u]:
            if u not in adj[v]:
                adj[v][u] = adj[u][v]
                G["m"] = G["m"] + 1
    return G

def checkCycle(G, node_list):
    cost = 0.0
    for i in range(len(node_list)):
        x = node_list[i]
        y = node_list[ (i + 1) % len(node_list) ]
        if y not in G["adj"][x]:
            return False, "Reported nodes do not form a cycle."
        cost +=  G["adj"][x][y]
    return True, cost

############################################################
# Functions that Make New Graphs
############################################################


def emptyGraph(n):
    G = {"n": n, "m":0, "adj": {}}
    for i in range(n):
        G["adj"][i] = {}
    return G

def cycleGraph(n):
    G = emptyGraph(n)
    for i in range(n):
        G["adj"][i][(i + 1) % n] = True
        G["adj"][i][(i - 1) % n] = True
        G["m"] = G["m"] + 2
    return G

def completeGraph(n, label = None):
    if label == None:
        default_val = True
    else:
        default_val = label
    G = emptyGraph(n)
    for i in range(n):
        for j in range(i+1,n):
            addDirEdge(G, i, j, label = default_val)
            addDirEdge(G, j, i, label = default_val)
    return G

def randomERGraph(n, p, seed = None):
    if seed == None:
        rng = np.random.default_rng()
    else:
        rng = np.random.default_rng(seed)
    G = emptyGraph(n)
    for i in range(n):
        for j in range(i+1,n):
            # This iterates over all pairs i,j where j>i
            if (rng.random() < p): #This happens with probability p (assuming good rng)
                G["adj"][i][j] = True
                G["adj"][j][i] = True
                G["m"] = G["m"] + 2
    return G

            
def sampleWoutReplace(n,d, exclude, rng):
    # generate d numbers without replacement from {0,...,n-1} - {exclude}.
    sample = [exclude]
    for j in range(d):
        nbr = rng.integers(n)
        while nbr in sample:
            nbr = rng.integers(n)
        sample.append(nbr)
    return sample[1:]


def dag(n,d):
    # creates a dag with roughly d out-degree
    # all edges go from smaller index to larger index
    rng = np.random.default_rng(None)
    G =  emptyGraph(n)
    for i in range(n-2):
        for k in range(d):
            j = rng.integers(low=i+1, high=n)
        addDirEdge(G,i,j)
    return G

def randomDigraphDegreeBound(n, d, seed=None):
    # Generate a random graph on n vertices where each vertex 
    # is given d random incoming and outgoing edges
    # (repeat edges are ignored)
    rng = np.random.default_rng(seed)
    G = emptyGraph(n)
    for i in G["adj"]:
        out_list = sampleWoutReplace(n,d,i, rng)
        for nbr in out_list:
            addDirEdge(G, i, nbr, label = rng.random())
        in_list = sampleWoutReplace(n,d,i, rng)
        for nbr in in_list:
            addDirEdge(G, nbr, i, label = rng.random())
    return G

def randomSignedDiGraph(n,d,q, seed = None):
    # Generates a graph with d random incoming and outgoing edges per node. Each each is given sign -1 with probability q and 1 with probability 1-q. 
    G = randomDigraphDegreeBound(n, d, seed = seed)
    rng = np.random.default_rng(seed)
    for u in G["adj"]:
        for v in G["adj"][u]:
            if rng.random() < q:
                G["adj"][u][v] = -1
            else:
                G["adj"][u][v] = +1
    return G

def oneNegCycle(n):
    G = emptyGraph(n)
    for i in range(n):
        addDirEdge(G, i , (i+1) % n, label = 1)
    G["adj"][0][1] = -n
    return G

############################################################
# Traversals
############################################################

def BFS(G, s):
    # G is a dictionary with keys "n", "m", "adj" representing an unweighted graph
    # G["adj"][u][v] is True if (u,v) is present. Otherwise, v is not in G["ad"][u].
    distances = {}
    finalized = {} # set of discovered nodes
    parents = {} # lists parent of node in SP tree
    layers = [[] for d in range(G["n"])] # lists of nodes at each distance.
    Q = queue.Queue()
    distances[s] = 0
    parents[s] = None
    Q.put(s)
    while not(Q.empty()): #Q not empty
        u = Q.get()
        if u not in finalized: #if u was already finalized, ignore it.
            finalized[u] = True
            layers[distances[u]].append(u) 
            for v in G["adj"][u]:
                # record v's distance and parent and add v to the queue if  
                # this is the first path to v,  
                if (v not in distances): # first path to v
                    distances[v] = distances[u] + 1
                    parents[v] = u
                    Q.put(v)
    return distances, parents, layers

def DFS(G):
    # G is a dictionary with keys "n", "m", "adj" representing an unweighted graph
    # Returns an array of discovered, containing timestamps of when each node is encountered
    color = {}
    discovered = {}
    finished = {}
    parent = {}
    for u in G["adj"]:
        color[u] = "white"
        parent[u] = None
    timestamp = [0] #This is a list whose only element is the current value of the time stamp. 

    def DFSVisit(u,  G, timestamp, color, discovered, finished):
        # Only the first argument ever changes
        color[u] = "gray"
        timestamp[0] = timestamp[0] + 1
        discovered[u] = timestamp[0]
        for v in G["adj"][u]:
            if color[v] == "white":
                    parent[v] = u
                    DFSVisit(v,  G, timestamp, color, discovered, finished)
        color[u] = "black"
        timestamp[0] = timestamp[0] + 1
        finished[u] = timestamp[0]
        return

    for u in G["adj"]:
        if color[u] == "white":
            DFSVisit(u, G, timestamp, color, discovered, finished)
    return discovered, finished, parent


def dijkstra(G, s):
    # G is a dictionary with keys "n", "m", "adj" representing an *weighted* graph
    # G["adj"][u][v] is the cost (length / weight) of edge (u,v)
    # This algorithms finds least-costs paths to all vertices
    # Returns an array of distances (path costs) and parents in the lightest-paths tree.
    # Assumes nonnegative path costs
    distances = {} # actual distances
    finalized = {} # set of discovered nodes
    parents = {} # lists parent of node in SP tree
    Q = [] # empty priority queue. Use heappush(Q, (priorit, val)) to add. Use heappop(Q) to remove.
    distances[s] = 0
    parents[s] = None
    heapq.heappush(Q, (distances[s], s))
    while len(Q) > 0: #Q not empty
        (d, u) = heapq.heappop(Q)
        if u not in finalized: #if u was already finalized, ignore it.
            finalized[u] = True
            for v in G["adj"][u]:
                new_length = distances[u] + G["adj"][u][v]
                # update v's distance (and parent and priority queue) if  
                # either this is the first path to v 
                # or we have found a better path to v
                if ((v not in distances) or (new_length < distances[v] )):
                    distances[v] = new_length
                    parents[v] = u
                    # add a copy of v to the queue with priority distances[v]
                    heapq.heappush(Q, (distances[v], v)) 
    return distances, parents