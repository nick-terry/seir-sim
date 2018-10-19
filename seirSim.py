# -*- coding: utf-8 -*-
"""
seir-sim

Simulates dynamics of infectious diseases on a contact network using a 
discretized SEIR model
"""

from random import random,choice
import numpy as np

class SeirSim():
    """
    Creates a simulation for the SEIR model for infectious disease
    
    Parameters:
        G: The contact network used to run the simulation (networkX graph)
        expRate(: The rate of conversion from Susceptible to Exposed
        infRate: The rate of conversion from Exposed to Infected
        recRate: The rate of conversion from Infected to Removed/Recovered
        tallyFuncs (optional): Functions to be used as tally statistics
            during the simulation run
        logSim (optional): Whether or not to record the simulations results
    """
    nodeStateDict = {0:'Susceptible',
                     1:'Exposed',
                     2:'Infectious',
                     3:'Recovered'}
    
    def __init__(self,G,expRate,infRate,recRate,
                 rngSeed=None,
                 tallyFuncs=None,
                 logSim=False):
        self.G = G
        self.expRate = expRate
        self.infRate = infRate
        self.recRate = recRate
        self.tallyFuncs = tallyFuncs
        self.logSim = logSim
        self.rngSeed = rngSeed
        
        self.nodeStates = np.zeros(len(G.nodes()))
    
        self.exposedList = []
    
        self.infectiousList = []
        self.siList = []
                    
        if not self.rngSeed is None:
            random.seed(self.rngSeed)
            
        #Mark a single node as infectious
        validSourceNodes = np.where(self.nodeStates==0)
        sourceNode = choice(validSourceNodes[0])
        self.sourceNode = sourceNode
        
        self.nodeStates[sourceNode] = 2
        self.infectiousList.append(sourceNode)
        for edge in G.edges(sourceNode):
            if self.nodeStates[edge[1]]==0:
                self.siList.append(edge)
        
        #A list of simulation state arrays
        if self.logSim:
            self.simStates = []
            self.simState = [self.nodeStates.copy(),self.siList.copy()]
            self.simStates.append(self.simState)
    
        self.t = 0
        
        self.useTally = not tallyFuncs is None
        
        if (self.useTally):
            self.tallyStats = []
            self.recordTallyStats()
    
    def __copy__(self):
        return SeirSim(self.G,self.expRate,self.infRate,self.recRate,
                 rngSeed=self.rngSeed,
                 tallyFuncs=self.tallyFuncs,
                 logSim=self.logSim)
    
    def transitionSE(self):
        '''
        S->E algorithm
        '''
        #Choose random edge (i,j) between SI
        edge = choice(self.siList)
        
        if (edge[0] in self.infectiousList):
            newExposedNode = edge[1]
        else:
            newExposedNode = edge[0]
        
        #Remove all edges from SI list which contain the newly exposed node
        for edge in self.siList.copy():
            if newExposedNode in edge:
                self.siList.remove(edge)
        
        #Add new exposed node to exposed list, mark as exposed in state array
        self.exposedList.append(newExposedNode)
        self.nodeStates[newExposedNode] = 1
        
    def transitionEI(self):
        '''
        E->I algorithm
        '''
        #Choose a random exposed node to become infectious
        newInfectedNode = choice(self.exposedList)
        
        #Remove node from exposed list and add to infected list
        self.exposedList.remove(newInfectedNode)
        self.infectiousList.append(newInfectedNode)
        
        #Check each edge (i,j) of the newly infected node i
        for edge in self.G.edges(newInfectedNode):
                #If j is infected, remove the edge from SI list

                if (self.nodeStates[edge[1]]==2):
                    if ((edge[1],newInfectedNode) in self.siList):    
                        self.siList.remove((edge[1],newInfectedNode))
                #If j is susceptible, add the edge to the SI list
                elif (self.nodeStates[edge[1]]==0):
                    self.siList.append((edge[1],newInfectedNode))
                    
        self.nodeStates[newInfectedNode] = 2
        
    def transitionIR(self):
        #I-R algorithm
        #Choose a random infectious node to remove
        newRemovedNode = choice(self.infectiousList)
        #Remove newly removed node from infectious list, update state array
        self.infectiousList.remove(newRemovedNode)
        self.nodeStates[newRemovedNode] = 3
        #Remove edges (i,j) from SI list, where i is the newly removed node
        for edge in self.siList.copy():
            if newRemovedNode in edge:
                self.siList.remove(edge)
        
    def simulate(self,policies=None):
        '''
        Runs the SEIR simulation
        
        Returns:
        simStates: An array which records the state of each node during
            each time step of the simulation
        tallyResults: A 2D list containing the tally results
        '''
        #print('Beginning simulation...')
        self.policies = policies
            
        #Execute any policies with runOnInit=True 
        if not self.policies is None:
            for policy in self.policies:
                if not policy is None and policy.runOnInit:
                    policy.execute(self)
        
        #Simulation tick loop
        while (len(self.siList) + len(self.exposedList) + len(self.infectiousList) > 0):
        
            self.t = self.t + 1 
        
            #True if an exposure event is possible
            expPossible = 1 if len(self.siList) > 0 else 0
            #True if an recovery event is possible
            recPossible = 1 if len(self.infectiousList) > 0 else 0
            #True if an infection event is possible
            infPossible = 1 if len(self.exposedList) > 0 else 0
                
            
            #Calculate transition probabilities
            #We exclude terms from the denominator if the corresponding event is
            #not possible
            probDenom = expPossible*self.expRate*len(self.siList)+infPossible*self.infRate*len(self.exposedList)+recPossible*self.recRate*len(self.infectiousList)
            probSE = expPossible*self.expRate*len(self.siList)/probDenom
            probEI = infPossible*self.infRate*len(self.exposedList)/probDenom
            
            x = random()
    
            #Produce transition event
            if (x < probSE):
                self.transitionSE()
            
            elif (x >= probSE and x < probEI+probSE):
                self.transitionEI()
            
            elif (recPossible):
                self.transitionIR()
            
            self.simState = [self.nodeStates.copy(),self.siList.copy()]
            
            if self.useTally:
                self.recordTallyStats()
                
            if self.logSim:
                self.simStates.append(self.simState)
            
            if not self.policies is None:
                for policy in self.policies:
                    if not policy is None and policy.runEachTimestep:
                        policy.execute(self)
            
        return self.simStates,np.array(self.tallyStats)


    def recordTallyStats(self):
        '''
        Helper method thats calls all tally functions passed to the simulation
        '''
        tResults = []
        for f in self.tallyFuncs:
            tResults.append(f(self.simState))
        self.tallyStats.append(tResults)
        return self.tallyStats