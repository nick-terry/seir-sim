# -*- coding: utf-8 -*-
"""
@author: Nick Terry

A test script for seir-sim that generates a random contact network and runs the
simulation
"""
import matplotlib.pyplot as plt
import numpy as np

from experiment import Experiment
from seirSim import SeirSim
from policies import vaccinateTopNDIL,vaccinateNRandom
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

DILPolicy = vaccinateTopNDIL(200)
randomVaccPolicy = vaccinateNRandom(200)

simulation = SeirSim(G,exposureRate,infectionRate,recoveryRate,
                     logSim=True, 
              tallyFuncs=[numNodesInState(0),numNodesInState(1),
                          numNodesInState(2),numNodesInState(3)])

exper = Experiment(simulation)

logList,statsList = exper.compare([[None],[randomVaccPolicy],[DILPolicy]])

def plotResults(axes,log,stats,title=None):
    for i in range(4):
        axes.plot(range(len(log)),stats[:,i])
    axes.legend(['Susceptible','Exposed','Infected','Removed'],loc='upper right')
    plt.xlabel('t')
    plt.ylabel('Population')
    if type(title) == str:
        plt.title(title)

def plotExperiment(experiment,titles):
    f = plt.figure(figsize=(10,3))
    axesList = []
    numSims = len(experiment.stats)
    for i in range(numSims):
        axesList.append(f.add_subplot('1{0}{1}'.format(numSims,i)))
        plotResults(axesList[i],experiment.logs[i],experiment.stats[i],
                    titles[i])

titles = ['No Vaccination','Random Vaccination','DIL-Ranked Vaccination']
plotExperiment(exper,titles)
