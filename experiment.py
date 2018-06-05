# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 19:09:01 2018

@author: pnter
"""
import copy

class Experiment():
    '''
    Experiments allow for the comparison of different policies on the same
    network/simulation. Random number generation will use the same seed in
    each simulation to allow for direct comparison of results
    '''
    def __init__(self,simulation):
        self.simulation = simulation
        self.stats = []
        self.logs = []
        
    def compare(self,policySets):
        '''
        Arguments:
            policySets: an array where each entry is an list of policies to use
                for a single simulation
        '''
        for policies in policySets:
            sim = copy.copy(self.simulation)
            logs,stats = sim.simulate(policies)
            self.logs.append(logs)
            self.stats.append(stats)
        
        return self.logs,self.stats