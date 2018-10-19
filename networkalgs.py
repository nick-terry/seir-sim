# -*- coding: utf-8 -*-
"""
Created on Thu May 24 21:46:25 2018

@author: pnter
"""
import numpy as np
import networkx as nx
from random import random,choice
from math import log1p
from scipy.stats import norm

def DIL(G):
    '''
    Computes DIL importance statistic for all nodes in graph G
    '''
    
    def countThreeCycles(edge):
        '''
        Counts the number of 3-cycles containing edge e_ij
        
        cycles is a map from node to the set of cycles containing node
        '''
        #Check previously found cycles for e_ij 
        start = edge[1]
        end = edge[0]
        
        count = 0
        
        for cycle in cycles[start]:
            if end in cycle:
                count += 1
        
        startEdges = G.edges(start)
        endEdges = G.edges(end)
        maxCycles = min(len(startEdges),len(endEdges))
        
        #e_ij can have no more cycles than degree(j)
        if count == maxCycles:
            return count
        
        for edge1 in startEdges:
            middle = edge1[1]
            for edge2 in G.edges(middle):
                if edge2[1] == end and {start,middle,end} not in cycles[start]:
                    count += 1
                    for node in start,middle,end:
                        if not {start,middle,end} in cycles[node]: cycles[node].append({start,middle,end})
                    break
            if count == maxCycles:
                break
            
        return count
    
    def I(edge):
        '''
        Computes the importance of an edge e_ij
        '''
        p = countThreeCycles(edge)
        u = (G.degree(edge[0])-p-1)*(G.degree(edge[1])-p-1)
        return u/(p/2+1)
    
    def W(i,j):
        '''
        Computes contribution of v_j to importance of v_i
        
        Assumes I_ij has already been computed
        '''
        degi = G.degree(i)
        degj = G.degree(j)
        if degi - degj == 0:
            return 0
        else:
            return edgeToIWArray[i,j,0]*(degi-1)/((degi+degj)-2)
    
    def L(i):
        return G.degree(i)+sum(edgeToIWArray[edge[0],edge[1],1] for edge in G.edges(i))
    
    cycles = [[] for i in G.nodes()]
    edgeToIWArray = np.zeros((len(G.nodes()),len(G.nodes()),2))

    for edge in G.edges():
        i,j = edge
        edgeToIWArray[i,j,0] = I(edge)
        edgeToIWArray[i,j,1] = W(i,j)
    
    nodeToL = np.zeros(len(G.nodes()))
    for node in G.nodes():
        nodeToL[node] = L(node)
    
    return nodeToL

def generateRandomGraph(numNodes, prob=None):
    '''
    Generates a random graph G with given number of Nodes and either the given
    probability probability p in the Erdos-Renyi G(n,p) model, or nodes with 
    the given degree distribution. 
    
    
    Parameters:
        numNodes: The number of nodes, n, on G
        prob (optional): The probability parameter p in the Erdos-Renyi random
            G(n,p) model. If prob is not provided, p=ln(n)/n will be used,
            as this is the threshold value for which G will almost surely be
            connected
    
    Returns:
        G: A networkx graph
    '''
    #print('Generating a random graph with {0} nodes...'.format(numNodes))
    
    if (prob is None):
        prob = log1p(numNodes)/numNodes
        
    G = nx.fast_gnp_random_graph(numNodes,prob)
    while not nx.is_connected(G):
        G = nx.fast_gnp_random_graph(numNodes,prob)
        
    return G

def sampleF(F):
    '''
    Generate a random realization from the cumulative distribution function F
    '''
    x = random()
    for i in range(F.size):
        if (F[i,1] > x):
            return int(F[i,0])    
        
def twoStepHeuristic(G,n1,n2):
    '''
    Use the two step heuristic algorithm described by K. Avrachenkov et al to
    search for the n nodes with highest degree.
    
    Perform the two step heuristic algorithm with
    n1 chosen for optimal performance according to the paper
    
    Parameters:
        G: The graph to perform the heuristic search on
        n1: The number of nodes to sample from G for the 1st step
        n2: The culling parameter for the 2nd step of the TSH algorithm
    
    Returns:
        topNodes: A list of the top-k nodes of G as found by the TSH algorithm
    '''
    nodes = G.nodes()
    S = np.zeros([len(G.nodes())])
    for i in range(n1):
        v = choice(nodes)
        for neighbor in G.neighbors(v):
            S[neighbor] += 1
    topNodes = np.argsort(S)[-n2:]

    return topNodes

def acquaintanceN(G,n):
    '''
    Use the acquaintance algorithm to select n nodes to vaccinate. This process
    biases the selected nodes to be of higher degree. Based on
    paper by Cohen, Havlin and ben Avraham.
    
    Parameters:
        G: The graph to perform the heuristic search on
        n: The number of nodes to sample
    
    Returns:
        topNodes: A list of nodes to vaccinate
    '''
    nodes = G.nodes()
    selectedNodes = []
    
    def surveyNodeForNeighbor(nodes):
        node = choice(nodes)
        neighbors = G.neighbors(node)
        neighbor = choice(neighbors)
        neighbors.remove(neighbor)
        
        while neighbor in selectedNodes and len(neighbors)>0:
            neighbor = choice(neighbors)
            neighbors.remove(neighbor)
        
        if neighbor in selectedNodes:
            nodes.remove(node)
            return surveyNodeForNeighbor(nodes)
        else:
            return(neighbor)
        
    for i in range(n):
        nodeToVacc = surveyNodeForNeighbor(nodes)
        selectedNodes.append(nodeToVacc)
            
    return selectedNodes
        
        