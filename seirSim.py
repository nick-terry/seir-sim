# -*- coding: utf-8 -*-
"""
seir-sim

Simulates dynamics of infectious diseases on a contact network using a 
discretized SEIR model
"""

import networkx as nx
from random import random,choice
import numpy as np

def generateRandomGraph(numNodes, numEdges):
    G = nx.Graph()

    for i in range(numNodes):
        G.add_node(i)


    for i in range(numEdges):
        u,v = choice(range(numNodes)),choice(range(numNodes))
        edges = G.edges()
        
        while u==v or (u,v) in edges or (v,u) in edges:
            u,v = choice(range(numNodes)),choice(range(numNodes))    
        G.add_edge(u,v)
    
    return G

def seirSim(G,expRate,infRate,recRate,tallyFuncs=None,logSim=False):
    """
    Runs a discretized SEIR simulation on a contact network
    
    Parameters:
        G: The contact network used to run the simulation (networkX graph)
        expRate(: The rate of conversion from Susceptible to Exposed
        infRate: The rate of conversion from Exposed to Infected
        recRate: The rate of conversion from Infected to Removed/Recovered
        tallyFuncs (float)(optional): Functions to be used as tally statistics
        during the simulation run
            
    Returns:
        tallyResults: A 2D list containing the tally results
    
    """
    #Create state array for nodes (0=Susceptible, 1=Exposed, 2=Infectious, 3=Recovered)
    nodeStates = np.zeros(len(G.nodes()))
    
    useTally = not tallyFuncs is None
    tallyStats = []
    
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
        simStates.append(nodeStates.copy())
    
    t = 0
    
    
    print('init::::::'+str(len(infectiousList)))
    #Simulation tick loop (while there are still nodes which can be infected)
    while (len(siList) > 0):
        
        t = t + 1 
        
        #Calculate transition probabilities
        probDenom = expRate*len(siList)+infRate*len(exposedList)+recRate*len(infectiousList)
        probSE = expRate*len(siList)/probDenom
        probEI = infRate*len(exposedList)/probDenom
        
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
                
            #Add new exposed node to exposed list, mark as exposed in state array
            exposedList.append(newExposedNode)
            nodeStates[newExposedNode] = 1
                
        
        elif (x >= probSE and x < probEI):
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
        
        elif (len(infectiousList) > 0):
            #I-R algorithm
            #Choose a random infectious node to remove
            newRemovedNode = choice(infectiousList)
            #Remove newly removed node from infectious list, update state array
            infectiousList.remove(newRemovedNode)
            nodeStates[newRemovedNode] = 3
            #Remove edges (i,j) from SI list, where i is the newly removed node
            for edge in siList:
                if newRemovedNode in edge:
                    siList.remove(edge)
        
        simState = [nodeStates,siList]
        if useTally:
            tallyStats = recordTallyStats(t,simState,tallyFuncs,tallyStats)
            
        if logSim:
            simStates.append(simState)
            
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

numNodes = 100
numEdges = 400

exposureRate = 10
infectionRate = 15
recoveryRate = .5

G = generateRandomGraph(numNodes, numEdges)

log,stats = seirSim(G,exposureRate,infectionRate,recoveryRate,logSim=True, 
              tallyFuncs=[numNodesInState(0),numNodesInState(1),
                          numNodesInState(2),numNodesInState(3)])

plt.plot(range(len(log)-1),stats[:,0],stats[:,1],stats[:,2],stats[:,3])
plt.legend(['Susceptible','Exposed','Infected','Recovered'],loc='upper left')
    
    
    