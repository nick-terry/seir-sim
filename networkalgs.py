# -*- coding: utf-8 -*-
"""
Created on Thu May 24 21:46:25 2018

@author: pnter
"""
import numpy as np
import networkx as nx
from random import random,choice

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
        return edgeToIWArray[i,j,0]*(degi-1)/(degi+G.degree(j)-2)
    
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

def generateRandomGraph(numNodes, numEdges=None, degreeDist=None):
    '''
    Generates a random graph G with given number of Nodes and either the given
    number of edges, or nodes with the given degree distribution
    
    Parameters:
        numNodes: The number of nodes on G
        numEdges (optional): The number of edges on G
        degreeDist (optional): The cumulative distribution function F giving the
            probability that a node N has degree D <= x. Takes the form of a 
            list of a 2D numpy array, where the 0th column are values x which 
            can be taken on by D, and 1st column are the P(D<=x)
    
    Returns:
        G: A networkx graph
    '''
    print('Generating a random graph with {0} nodes...'.format(numNodes))
    G = nx.Graph()

    for i in range(numNodes):
        G.add_node(i)
        
    if (numEdges is None and degreeDist is None):
        raise Exception('You must provide a number or edges or degree distribution')
        
    elif (not numEdges is None):
        #Choose 2 nodes randomly from a uniform distribution and create an
        #edge between them
        for i in range(numEdges):
            u,v = choice(range(numNodes)),choice(range(numNodes))
            edges = G.edges()
            
            while u==v or (u,v) in edges or (v,u) in edges:
                u,v = choice(range(numNodes)),choice(range(numNodes))    
            G.add_edge(u,v)
            
    if (not degreeDist is None):
        #Algorithm for generating a random graph with N nodes, each with degree
        #belonging to distribution D
        nodes = G.nodes().copy()
        
        #Generate the degree of each node from D
        degrees = []
        for ind in range(len(nodes)):
            degrees.append(sampleF(degreeDist))
            
        #A 2D array where array[node] = [target degree, degree of the node]
        nodeToEdgesNeeded = np.array([degrees,  degrees]).transpose()

        #Check each node, and add edges until it has the correct degree
        while (len(nodes) > 0):
            #Choose a random node (uniformly)
            i = nodes[0]
            
            #Find all nodes j such that:
            # 1) (i,j) is not an edge in G
            # 2) degree(j) < x_j, where x_j is the target degree of j
            # 3) j != i
            inEdgesFn = np.vectorize(lambda x: x not in np.ravel(G.edges(i)))
            print(nodeToEdgesNeeded[:,0] != i)
            candidateJ = list(np.where(np.logical_and(nodeToEdgesNeeded[:,1] > 0, inEdgesFn(nodeToEdgesNeeded[:,0])))[0])
            if i in candidateJ:
                candidateJ.remove(i)  
            
            while nodeToEdgesNeeded[i,1] > 0:
                #Randomly choose a candidate node j to create a new edge (i,j)
                j = int(choice(candidateJ))
                G.add_edge(i,j)
                #Update remaining number of edges needed by nodes i and j
                nodeToEdgesNeeded[i,1] = nodeToEdgesNeeded[i,1] - 1
                nodeToEdgesNeeded[j,1] = nodeToEdgesNeeded[j,1] - 1
                
            nodes.remove(i)
        
    return G

def sampleF(F):
    '''
    Generate a random realization from the cumulative distribution function F
    '''
    x = random()
    for i in range(F.size):
        if (F[i,1] > x):
            return int(F[i,0])    