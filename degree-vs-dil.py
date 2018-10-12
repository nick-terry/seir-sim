# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 22:32:32 2018

@author: pnter
"""

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
from plotting import plotExperiment


numNodes = 10000

G = generateRandomGraph(numNodes)

#Use parameters estimated from Ebola epidemics in Sierra Leone
#Source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4169395/
exposureRate = .45 #Beta
infectionRate = 1/5.3 #Sigma
recoveryRate = 1/5.61 #Gamma

numVaccinated = 2000

DILPolicy = vaccinateTopNDIL(numVaccinated)
DegPolicy = vaccinateTopNDegree(numVaccinated)

simulation = SeirSim(G,exposureRate,infectionRate,recoveryRate,
                     logSim=True, 
              tallyFuncs=[numNodesInState(0),numNodesInState(1),
                          numNodesInState(2),numNodesInState(3)])

exper = Experiment(simulation)

logList,statsList = exper.compare([[DegPolicy],[DILPolicy],])

titles = ['Degree-Ranked Vaccination','DIL-Ranked Vaccination']

plotExperiment(exper,titles)