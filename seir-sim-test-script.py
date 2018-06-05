# -*- coding: utf-8 -*-
"""
@author: Nick Terry

A test script for seir-sim that generates a random contact network and runs the
simulation
"""
import matplotlib.pyplot as plt
import numpy as np

from seirSim import SeirSim
from policies import vaccinateTopNDIL
from tallyFuncs import numNodesInState
from networkalgs import generateRandomGraph


numNodes = 1000
degreeDist = np.array([[1, .1],
              [2,.4],
              [3,.6],
              [4,.7],
              [5,1]])

G = generateRandomGraph(numNodes, degreeDist=degreeDist)

exposureRate = 10
infectionRate = 3
recoveryRate = .5

vaccPolicy = vaccinateTopNDIL(200)

simulation = SeirSim(G,exposureRate,infectionRate,recoveryRate,
                     policiesList=[vaccPolicy],
                     logSim=True, 
              tallyFuncs=[numNodesInState(0),numNodesInState(1),
                          numNodesInState(2),numNodesInState(3)])

log,stats = simulation.simulate()

for i in range(4):
    plt.plot(range(len(log)),stats[:,i])
plt.legend(['Susceptible','Exposed','Infected','Removed'],loc='upper right')
plt.xlabel('t')
plt.ylabel('Population')