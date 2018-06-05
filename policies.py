# -*- coding: utf-8 -*-
"""
Contains policies for use with seir-sim

Policies are inputs to the simulation which affects its evolution
"""
import numpy as np
from random import choice

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
            '''
            Ensure source node is not vaccinated. This ensures the simulation
            can actually run, and is justified by the fact that giving a
            vaccination to someone already infected by a disease does not cure
            the disease
            '''
            nodeDIL[simulation.sourceNode] = 0
            
            for i in range(self.N):
                topNode = np.argmax(nodeDIL)
                nodeDIL[topNode] = 0
                simulation.nodeStates[topNode] = 3
                
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
            '''
            Ensure source node is not vaccinated. This ensures the simulation
            can actually run, and is justified by the fact that giving a
            vaccination to someone already infected by a disease does not cure
            the disease
            '''
            nodes.remove(simulation.sourceNode)
            
            for i in range(self.N):
                node = choice(nodes)
                simulation.nodeStates[node] = 3
                nodes.remove(node)