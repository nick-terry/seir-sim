# -*- coding: utf-8 -*-
"""
seir-sim

Simulates dynamics of infectious diseases on a contact network using a 
discretized SEIR model
"""

import networkx as nx
from random import random,choice
import numpy as np

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
            # 2) degree(j) < x_j, where x_j is a realization of D
            
            inEdgesFn = np.vectorize(lambda x: x not in np.ravel(G.edges(i)))
            candidateJ = np.where(np.logical_and(nodeToEdgesNeeded[:,1] > 0, inEdgesFn(nodeToEdgesNeeded[:,0])))[0]
            
            while nodeToEdgesNeeded[i,1] > 0:
                #Randomly choose a candidate node j to create a new edge (i,j)
                j = choice(candidateJ)
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
            return F[i,0]


def seirSim(G,expRate,infRate,recRate,tallyFuncs=None,logSim=False):
    """
    Runs a discretized SEIR simulation on a contact network
    
    Parameters:
        G: The contact network used to run the simulation (networkX graph)
        expRate(: The rate of conversion from Susceptible to Exposed
        infRate: The rate of conversion from Exposed to Infected
        recRate: The rate of conversion from Infected to Removed/Recovered
        tallyFuncs (optional): Functions to be used as tally statistics
            during the simulation run
        logSim (optional): Whether or not to record the simulations results
            
    Returns:
        tallyResults: A 2D list containing the tally results
    
    """
    #Create state array for nodes (0=Susceptible, 1=Exposed, 2=Infectious, 3=Recovered)
    nodeStates = np.zeros(len(G.nodes()))
    
    exposedList = []
    
    infectiousList = []
    siList = []
    
    #Mark a single node as infectious
    sourceNode = choice(G.nodes())
    nodeStates[sourceNode] = 2
    infectiousList.append(sourceNode)
    for edge in G.edges(sourceNode):
        siList.append(edge)
    
    #A list of simulation state arrays
    if logSim:
        simStates = []
        simState = [nodeStates.copy(),siList.copy()]
        simStates.append(simState)
    
    t = 0
    
    useTally = not tallyFuncs is None
    
    if (useTally):
        tallyStats = []
        tallyStats = recordTallyStats(t,simState,tallyFuncs,tallyStats)
    
    print('Beginning simulation...')
    #Simulation tick loop
    while (len(siList) + len(exposedList) + len(infectiousList) > 0):
        
        t = t + 1 
        
        #True if an exposure event is possible
        expPossible = 1 if len(siList) > 0 else 0
        #True if an recovery event is possible
        recPossible = 1 if len(infectiousList) > 0 else 0
        #True if an infection event is possible
        infPossible = 1 if len(exposedList) > 0 else 0
            
        
        #Calculate transition probabilities
        #We exclude terms from the denominator if the corresponding event is
        #not possible
        probDenom = expPossible*expRate*len(siList)+infPossible*infRate*len(exposedList)+recPossible*recRate*len(infectiousList)
        probSE = expPossible*expRate*len(siList)/probDenom
        probEI = infPossible*infRate*len(exposedList)/probDenom
        
        x = random()
    
        #Produce transition event
        if (x < probSE):
            #S->E algorithm
            #Choose random edge (i,j) between SI
            edge = choice(siList)
            
            if (edge[0] in infectiousList):
                newExposedNode = edge[1]
            else:
                newExposedNode = edge[0]
            
            #Remove all edges from SI list which contain the newly exposed node
            for edge in siList.copy():
                if newExposedNode in edge:
                    siList.remove(edge)
            
            #Add new exposed node to exposed list, mark as exposed in state array
            exposedList.append(newExposedNode)
            nodeStates[newExposedNode] = 1

        
        elif (x >= probSE and x < probEI+probSE):
            #E->I algorithm
            #Choose a random exposed node to become infectious
            newInfectedNode = choice(exposedList)
            
            #Remove node from exposed list and add to infected list
            exposedList.remove(newInfectedNode)
            infectiousList.append(newInfectedNode)
            
            #Check each edge (i,j) of the newly infected node i
            for edge in G.edges(newInfectedNode):
                    #If j is infected, remove the edge from SI list

                    if (nodeStates[edge[1]]==2):
                        if ((edge[1],newInfectedNode) in siList):    
                            siList.remove((edge[1],newInfectedNode))
                    #If j is susceptible, add the edge to the SI list
                    elif (nodeStates[edge[1]]==0):
                        siList.append((edge[1],newInfectedNode))
                        
            nodeStates[newInfectedNode] = 2
        
        elif (recPossible):
            #I-R algorithm
            #Choose a random infectious node to remove
            newRemovedNode = choice(infectiousList)
            #Remove newly removed node from infectious list, update state array
            infectiousList.remove(newRemovedNode)
            nodeStates[newRemovedNode] = 3
            #Remove edges (i,j) from SI list, where i is the newly removed node
            for edge in siList.copy():
                if newRemovedNode in edge:
                    siList.remove(edge)
        
        simState = [nodeStates.copy(),siList.copy()]
        
        if useTally:
            tallyStats = recordTallyStats(t,simState,tallyFuncs,tallyStats)
            
        if logSim:
            simStates.append(simState)
        '''
        print('Time: '+str(t))
        print(exposedList)
        '''
    return simStates,np.array(tallyStats)


def recordTallyStats(t,simState,fList,results):
    '''
    Helper function thats calls all tally functions passed to the simulation
    
    Parameters:
        t: time index of the simulation
        simState: state of the simulation at time t
        fList: list of functions which take a simState as input and return a 
            scalar
        results: the list of tally statistics maintained by the simulation
    '''
    tResults = []
    for f in fList:
        tResults.append(f(simState))
    results.append(tResults)
    return results
                            
def numNodesInState(state):
    '''
    Decorator function for creating tally statistics for nodes in
    a certain state
    '''
    def numNodesInState_(simState):
        '''
        Tally statistic which computes the number of nodes in state at time t
        '''
        n = len(simState[0])
        return np.sum(np.where(simState[0]==state,np.ones(n),np.zeros(n)))
    
    return numNodesInState_