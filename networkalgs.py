# -*- coding: utf-8 -*-
"""
Created on Thu May 24 21:46:25 2018

@author: pnter
"""
import networkx
import numpy as np

def DIL(G):
    
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