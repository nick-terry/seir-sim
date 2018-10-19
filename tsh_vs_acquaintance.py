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
from experiment import Experiment
from seirSim import SeirSim
from policies import vaccinateTopNTSH,vaccinateNAcquaintance,vaccinateNRandom
from tallyFuncs import numNodesInState
from networkalgs import generateRandomGraph
from plotting import plotExperiment
import numpy as np
import matplotlib.pyplot as pl
from scipy.stats import norm

numNodes = 1000

#Use parameters estimated from Ebola epidemics in Sierra Leone
#Source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4169395/
exposureRate = .45 #Beta
infectionRate = 1/5.3 #Sigma
recoveryRate = 1/5.61 #Gamma

numVaccinated = 500
numIterations = 100

RandomPolicy = vaccinateNRandom(numVaccinated)
TSHPolicy = vaccinateTopNTSH(numVaccinated)
AcqPolicy = vaccinateNAcquaintance(numVaccinated)
policies = [[RandomPolicy],[TSHPolicy],[AcqPolicy]]
#Store statistics on proportion of population infected (Average,Var)
proportionInfected = np.zeros([numIterations,len(policies),2])
#Store count of # nodes infected in each simulation
numInfected = np.zeros([numIterations,len(policies)])

for iteration in range(numIterations):
    print("Beginning Simulation {0}".format(iteration))
    G = generateRandomGraph(numNodes)
    
    simulation = SeirSim(G,exposureRate,infectionRate,recoveryRate,
                         logSim=True, 
                  tallyFuncs=[numNodesInState(0),numNodesInState(1),
                              numNodesInState(2),numNodesInState(3)])
    
    exper = Experiment(simulation)
    
    logList,statsList = exper.compare(policies)
    abort = False
    for policy in range(len(policies)):
        propInfectedMean = np.mean(statsList[policy][:,2]/numNodes)
        propInfectedVar = np.var(statsList[policy][:,2]/numNodes)
        proportionInfected[iteration,policy,:] = np.array([propInfectedMean,propInfectedVar])
        
        numRemoved = np.sum(np.where(logList[policy][-1][0]==3,
                                   np.ones(numNodes),np.zeros(numNodes)))
        numInfected[iteration,policy] = numRemoved - numVaccinated

titles = ['Random Immunization','Two-Step Heuristic Immunization','Acquaintance Immunization']

#plotExperiment(exper,titles)
ax = pl.subplot(1,1,1)
ax.plot(np.array(range(numIterations)),proportionInfected[:,0,0],' . ',
        color='blue',linestyle='None')
ax.plot(np.array(range(numIterations)),proportionInfected[:,1,0],' . ',
        color='orange',linestyle='None')
ax.plot(np.array(range(numIterations)),proportionInfected[:,2,0],' . ',
        color='green',linestyle='None')

box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])
    
# Put a legend below current axis
ax.legend(titles,loc='upper center', bbox_to_anchor=(0.5, -0.15),
          fancybox=True, shadow=True, ncol=len(policies))

pl.ylabel('Mean Proportion of Population Infected')
pl.xlabel('Simulation Replication #')

pl.title('Comparison of Random, TSH, and Acquaintance Immunization,\n K={0}'.format(
        numVaccinated))

pl.show()




#Compute mean, variance, 90% CI for each policy
for i in range(len(policies)):
    mean = np.mean(proportionInfected[:,i,0])
    var= np.var(proportionInfected[:,i,0])
    ciHalfwidth = norm.ppf(.9)*np.sqrt(var/numIterations)
    statList = [mean,var,ciHalfwidth]
    roundFn = lambda x:round(x,3)
    printList = statList
    printList.insert(0,titles[i])
    print('{0}:\nMean:{1}, Variance:{2}, CI Halfwidth:{3}'.format(
            *printList))
    print('# Infected:\nMean: {0}, Variance: {1}\n'.format(
            np.mean(numInfected[:,i]),np.var(numInfected[:,i])))
    

    