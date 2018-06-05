# -*- coding: utf-8 -*-
"""
Contains policies for use with seir-sim

Policies are inputs to the simulation which affects its evolution
"""
import numpy as np

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
            topN = []
            
            for i in range(self.N):
                topNode = np.argmax(nodeDIL)
                topN.append(topNode)
                nodeDIL[topNode] = 0
            
            for node in topN:
                simulation.nodeStates[node] = 3