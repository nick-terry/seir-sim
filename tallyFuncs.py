# -*- coding: utf-8 -*-
import numpy as np

"""
Contains tally functions for computing statistics in seir-sim
"""
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