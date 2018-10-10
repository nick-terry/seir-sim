# -*- coding: utf-8 -*-
"""
@author: Nick Terry

A test script for seir-sim that generates a random contact network and runs the
simulation
"""
import matplotlib.pyplot as plt

from experiment import Experiment
from seirSim import SeirSim
from policies import vaccinateTopNDIL,vaccinateNRandom,vaccinateTopNDegree,vaccinateTopNTSH, vaccinateNAcquaintance
from tallyFuncs import numNodesInState
from networkalgs import generateRandomGraph
import networkViz as nv
from math import ceil,floor


numNodes = 1000

G = generateRandomGraph(numNodes)

exposureRate = 10
infectionRate = 3
recoveryRate = .5

numVaccinated = 200

DILPolicy = vaccinateTopNDIL(numVaccinated)
randomVaccPolicy = vaccinateNRandom(numVaccinated)
DegPolicy = vaccinateTopNDegree(numVaccinated)
tshPolicy = vaccinateTopNTSH(numVaccinated)
acquPolicy = vaccinateNAcquaintance(numVaccinated)

simulation = SeirSim(G,exposureRate,infectionRate,recoveryRate,
                     logSim=True, 
              tallyFuncs=[numNodesInState(0),numNodesInState(1),
                          numNodesInState(2),numNodesInState(3)])

exper = Experiment(simulation)

logList,statsList = exper.compare([[None],[randomVaccPolicy],
                                   [DegPolicy],[DILPolicy],
                                   [tshPolicy],[acquPolicy]])

def plotResults(axes,log,stats,title=None,legend=False):
    for i in range(4):
        axes.plot(range(len(log)),stats[:,i])
    if legend:
        axes.legend(['Susceptible','Exposed','Infected','Removed'],loc='upper right')
    plt.xlabel('t')
    plt.ylabel('Population')
    if type(title) == str:
        plt.title(title)

def plotExperiment(experiment,titles):
    f = plt.figure(figsize=(12,3))
    axesList = []
    numSims = len(experiment.stats)
    for i in range(numSims):
        axesList.append(f.add_subplot('{0}{1}{2}'.format(ceil(numSims/4),
                                      4,i)))
        plotResults(axesList[i],experiment.logs[i],experiment.stats[i],
                    titles[i])

titles = ['No Vaccination','Random Vaccination',
          'Degree-Ranked Vaccination','DIL-Ranked Vaccination',
          'TSH-Ranked Vaccination','Acquaintance Vaccination']
plotExperiment(exper,titles)