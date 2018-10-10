# -*- coding: utf-8 -*-
"""
Contains policies for use with seir-sim

Policies are inputs to the simulation which affects its evolution
"""
import numpy as np
from random import choice
from math import floor

import networkalgs

class Policy():
    '''
    An input to the simulation
    
    Parameters:
        runOnInit: True if the policy is applied on simulation initialization;
            False by default
        runEachTimestep: True if the policy is applied on each timestep of the
            simulation; False by default
    '''
    def __init__(self,runOnInit=False,runEachTimestep=False):
        self.runOnInit = runOnInit
        self.runEachTimestep = runEachTimestep
        
    def execute(self,simulation):
        return
    
class vaccinateTopNDIL(Policy):
        '''
        Vaccinate the top N nodes ranked by the DIL metric as defined in
        "Evaluating the importance of nodes in complex networks", Liu et al., 2016
        
        Vaccinating a node changes its state from Susceptible to Removed
        
        Parameters:
            N: The number of top nodes to vaccinate
        '''
        def __init__(self,N):
            Policy.__init__(self,runOnInit=True,runEachTimestep=False)
            self.N = N
            
        def execute(self,simulation):
            nodeDIL = networkalgs.DIL(simulation.G)
            
            for i in range(self.N):
                topNode = np.argmax(nodeDIL)
                nodeDIL[topNode] = 0
                simulation.nodeStates[topNode] = 3
                
            '''
            Ensure source node is not vaccinated. This ensures the simulation
            can actually run, and is justified by the fact that giving a
            vaccination to someone already infected by a disease does not cure
            the disease
            '''
            simulation.nodeStates[simulation.sourceNode] = 2
                
class vaccinateTopNDegree(Policy):
        '''
        Vaccinate the top N nodes ranked by degree
        
        Vaccinating a node changes its state from Susceptible to Removed
        
        Parameters:
            N: The number of top nodes to vaccinate
        '''
        def __init__(self,N):
            Policy.__init__(self,runOnInit=True,runEachTimestep=False)
            self.N = N
            
        def execute(self,simulation):
            nodeDeg = list(map(lambda x:simulation.G.degree(x),simulation.G.nodes()))
            
            for i in range(self.N):
                topNode = np.argmax(nodeDeg)
                nodeDeg[topNode] = 0
                simulation.nodeStates[topNode] = 3
                
            '''
            Ensure source node is not vaccinated. This ensures the simulation
            can actually run, and is justified by the fact that giving a
            vaccination to someone already infected by a disease does not cure
            the disease
            '''
            simulation.nodeStates[simulation.sourceNode] = 2
                
class vaccinateNRandom(Policy):
        '''
        Vaccinate N nodes at random
        
        Vaccinating a node changes its state from Susceptible to Removed
        
        Parameters:
            N: The number of top nodes to vaccinate
        '''
        def __init__(self,N):
            Policy.__init__(self,runOnInit=True,runEachTimestep=False)
            self.N = N
            
        def execute(self,simulation):
            nodes = list(range(len(simulation.G.nodes())))
            
            for i in range(self.N):
                node = choice(nodes)
                simulation.nodeStates[node] = 3
                nodes.remove(node)
                
            '''
            Ensure source node is not vaccinated. This ensures the simulation
            can actually run, and is justified by the fact that giving a
            vaccination to someone already infected by a disease does not cure
            the disease
            '''
            simulation.nodeStates[simulation.sourceNode] = 2
                
class vaccinateTopNTSH(Policy):
        '''
        Vaccinate the top N nodes ranked by the Two Step Heuristic algorithm
        
        Vaccinating a node changes its state from Susceptible to Removed
        
        Parameters:
            N: The number of top nodes to vaccinate
        '''
        def __init__(self,N):
            Policy.__init__(self,runOnInit=True,runEachTimestep=False)
            self.N = N
            
        def execute(self,simulation):
            topInd = networkalgs.twoStepHeuristic(simulation.G,self.N*2,self.N)
 
            for node in topInd:
                simulation.nodeStates[node] = 3
                
            '''
            Ensure source node is not vaccinated. This ensures the simulation
            can actually run, and is justified by the fact that giving a
            vaccination to someone already infected by a disease does not cure
            the disease
            '''
            simulation.nodeStates[simulation.sourceNode] = 2
            
class vaccinateNAcquaintance(Policy):
        '''
        Vaccinate N nodes chosen by the Acquaintance algorithm
        
        Vaccinating a node changes its state from Susceptible to Removed
        
        Parameters:
            N: The number of top nodes to vaccinate
        '''
        def __init__(self,N):
            Policy.__init__(self,runOnInit=True,runEachTimestep=False)
            self.N = N
            
        def execute(self,simulation):
            topInd = networkalgs.acquaintanceN(simulation.G,self.N)
            
            for node in topInd:
                simulation.nodeStates[node] = 3
                
            '''
            Ensure source node is not vaccinated. This ensures the simulation
            can actually run, and is justified by the fact that giving a
            vaccination to someone already infected by a disease does not cure
            the disease
            '''
            simulation.nodeStates[simulation.sourceNode] = 2